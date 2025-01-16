from firedrake import *

comm = COMM_WORLD
mesh = UnitDiskMesh(3)

V = FunctionSpace(mesh, 'CG', 1)

def base(mesh, k):
    x, y = SpatialCoordinate(mesh)
    theta = atan2(x, y)
    cos_theta = cos(Constant(k)*theta)
    return cos_theta # + 1j*sin_theta

n = 10
M_np = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        M_np[i, j] = assemble(base(mesh, i)*base(mesh, j)*ds(domain=mesh))

M_inv_np = np.linalg.inv(M_np)

D = PETSc.Mat().createDense(size=(n, n), comm=comm)
PETSc.Sys.syncPrint(f"[{comm.rank}/{comm.size}] D size: {D.getSizes()}")
PETSc.Sys.syncFlush()

coeff = np.zeros(n)
for i in range(n):
    for j in range(n):
        D.setValue(i, j, coeff[i]*M_inv_np[i, j])
D.assemble()


v = TestFunction(V)
bs = []
for i in range(n):
    bs.append(
        assemble(base(mesh, i)*v*ds)
    )

bc = DirichletBC(V, 0, 'on_boundary')
bc_nodes = bc.nodes[bc.nodes < V.dof_dset.sizes[1]]
ncl = bc_nodes.size
B = PETSc.Mat().createDense(size=((None, n), (ncl, None)), comm=comm)
rend = comm.scan(ncl)
rstart = rend - ncl
bc_nodes_only_index = np.arange(rstart, rend, dtype=np.int32)
for i in range(n):
    B.setValues(i, bc_nodes_only_index, bs[i].dat.data[bc_nodes])
B.assemble()
PETSc.Sys.syncPrint(f"[{comm.rank}/{comm.size}] B size: {B.getSizes()}")
PETSc.Sys.syncFlush()

C = D.ptap(B) # C = B^T D B
PETSc.Sys.syncPrint(f"[{comm.rank}/{comm.size}] C size: {C.getSizes()}")
PETSc.Sys.syncFlush()

n_local = V.dof_dset.sizes[1]
lgmap = V.local_to_global_map(bcs=None)
A = PETSc.Mat().createAIJ(size=((n_local, None), (n_local, None)), comm=comm)
A.setLGMap(lgmap)
gindex = lgmap.apply(bc_nodes)
gcols = comm.allgather((rstart, rend, gindex))
for i_rstart, i_rend, i_global in gcols:
    A.setValues(gindex, i_global, C.getValues(bc_nodes_only_index, range(i_rstart, i_rend))) 
A.assemble()
A.view()

PETSc.Sys.syncPrint(f"[{comm.rank}/{comm.size}] A size: {A.getSizes()}")
PETSc.Sys.syncFlush()