from firedrake import *
from firedrake.output import VTKFile as File

mu = 0.01
T = 5

N_S = 64
N_T = 512*T

tau = T/N_T
h = 1/N_S

mesh = RectangleMesh(N_S, N_S, 1, 1)

u_0 = as_vector((0, 0))
u_D = as_vector((1, 0))
f = as_vector([0, 0])

cell = mesh.ufl_cell()
tdim = cell.topological_dimension()

# Mini element: P1 X P1b
P1 = FiniteElement("CG", cell, 1)
B = FiniteElement("B", cell, tdim+1)
P1b = P1 + B # or P1b = NodalEnrichedElement(P1, B)

V_u = VectorFunctionSpace(mesh, P1b)
V_p = FunctionSpace(mesh, "CG", 1)
V = MixedFunctionSpace([V_u, V_p])

w = Function(V) # u and p
u, p = split(w)

v, q = TestFunctions(V)

w_nm1 = Function(V)
u_nm1, p_nm1 = w_nm1.subfunctions
u_nm1.rename('u_h') # for visualization in paraview
p_nm1.rename('p_h')

F = (
      Constant(1/tau)*inner(u - u_nm1, v)*dx
    + mu*inner(grad(u), grad(v))*dx
    + 1/2*inner(dot(grad(u), u/2), v)*dx
    - 1/2*inner(dot(grad(v), u/2), u)*dx
    - p*div(v)*dx
    + div(u)*q*dx
    - inner(f, v)*dx(domain=mesh)
)

bc1 = DirichletBC(V.sub(0).sub(0), 0, (1, 2))
bc2 = DirichletBC(V.sub(0).sub(1), 0, (3, 4))
bc3 = DirichletBC(V.sub(0).sub(0), 1, 4)  # upper boundary

nullspace = MixedVectorSpaceBasis(V, [V.sub(0), VectorSpaceBasis(constant=True, comm=mesh.comm)])

problem = NonlinearVariationalProblem(F, w, bcs=[bc1, bc2, bc3])  # F = 0
solver = NonlinearVariationalSolver(problem,
                                    options_prefix='ns',
                                    solver_parameters=None, # {'snes_converged_reason': None, 'snes_max_it': 100},
                                    nullspace=nullspace
                                   )

u_, p_ = w.subfunctions

output = VTKFile('pvd/ns-equation.pvd')

u_nm1.assign(0)
output.write(u_nm1, p_nm1, time=0)

# for i in range(N_T):
for i in ProgressBar('Time').iter(range(N_T)):
    t = tau*(i+1)
    solver.solve()
    w_nm1.assign(w)
    if (i+1)%32 == 0:
        output.write(u_nm1, p_nm1, time=t)


