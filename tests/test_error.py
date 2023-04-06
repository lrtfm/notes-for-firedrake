from firedrake import *
from firedrake.petsc import PETSc
from firedrake.solving_utils import KSPReasons

def _make_reasons(reasons):
    return dict([(getattr(reasons, r), r)
                 for r in dir(reasons) if not r.startswith('_')])

PCFailedReason = _make_reasons(PETSc.PC.FailedReason())

def printf(*args, **kwargs):
    PETSc.Sys.Print(*args, **kwargs)

def get_ksp_reason(solver):
    r = solver.snes.getKSP().getConvergedReason()
    return KSPReasons[r]

def get_pc_failed_reason(solver):
    pc = solver.snes.getKSP().getPC()
    r = pc.getFailedReason()
    return PCFailedReason[r]
    
rank, size = COMM_WORLD.rank, COMM_WORLD.size

opts = PETSc.Options()
N = opts.getInt('N', 32*size)

init_mesh = RectangleMesh(nx=N, ny=N, Lx=1, Ly=1)
meshes = MeshHierarchy(init_mesh, 1)
test_mesh = meshes[-1]
x, y = SpatialCoordinate(test_mesh)
f = sin(pi*x)*sin(pi*y)

V = FunctionSpace(test_mesh, 'CG', degree=1)
u, v = TrialFunction(V), TestFunction(V)
a = inner(grad(u), grad(v))*dx - inner(f, v)*dx
bc = DirichletBC(V, 0, sub_domain='on_boundary')

u_h = Function(V, name='u_h')
problem = LinearVariationalProblem(lhs(a), rhs(a), u_h, bcs=bc)

solver_parameters = {'ksp_type': 'cg',
                     'ksp_max_it': 4, 
                     'ksp_converged_reason': None, 
                     'ksp_error_if_not_converged': None,
                     'pc_type': 'none'}
solver = LinearVariationalSolver(problem, solver_parameters=solver_parameters, options_prefix='')

for i in range(3):
    printf(f"Loop i = {i}")
    try:
        solver.solve()
    except ConvergenceError:
        printf(f"  Error from Firedrake: solver did not converged: {get_ksp_reason(solver)}")
        printf(f"                                             PC : {get_pc_failed_reason(solver)}")

    except PETSc.Error as e:
        if e.ierr == 91:
            printf(f"  Error from PETSc: solver did not converged: {get_ksp_reason(solver)}")
            printf(f"                                         PC : {get_pc_failed_reason(solver)}")
        # We should terminate the process as suggested by Matt https://lists.mcs.anl.gov/pipermail/petsc-users/2023-March/048146.html
        raise
