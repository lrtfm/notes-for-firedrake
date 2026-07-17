from firedrake import *
from firedrake.petsc import PETSc

opts = PETSc.Options()
N = opts.getInt('N', default=8)
test_mesh = RectangleMesh(nx=N, ny=N, Lx=1, Ly=1)

x, y = SpatialCoordinate(test_mesh)

# Manufactured solution: u = x^2 - 1/3, f = -2,
# g_N = 2 on the right boundary (tag 2) and 0 elsewhere
u_exact = x**2 - 1/3
f = Constant(-2)
g_N = Constant(2)

V = FunctionSpace(test_mesh, 'CG', degree=1)
R = FunctionSpace(test_mesh, 'R', 0)
W = V * R
PETSc.Sys.Print(f'Number of Dofs: {W.dim()}')

u, mu = TrialFunctions(W)
v, eta = TestFunctions(W)

a = inner(grad(u), grad(v))*dx + inner(mu, v)*dx + inner(u, eta)*dx
L = inner(f, v)*dx + inner(g_N, v)*ds(2)

w_h = Function(W)
solve(a == L, w_h, options_prefix='test')

u_h, mu_h = w_h.subfunctions
PETSc.Sys.Print(f'L2 error: {errornorm(u_exact, u_h)}')

filename = 'pvd/u_h_neumann.pvd'
PETSc.Sys.Print(f'Write pvd file: {filename}')
VTKFile(filename).write(u_h)
