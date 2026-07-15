from firedrake import *

N = 8
test_mesh = RectangleMesh(nx=N, ny=N, Lx=1, Ly=1)
x, y = SpatialCoordinate(test_mesh)
u_exact = sin(pi*x)*sin(pi*y)
f = 2*pi**2*u_exact

V = FunctionSpace(test_mesh, 'CG', degree=1)
u, v = TrialFunction(V), TestFunction(V)
a = inner(grad(u), grad(v))*dx
L = inner(f, v)*dx
bc = DirichletBC(V, 0, 'on_boundary')

u_h = Function(V, name='u_h')
solve(a == L, u_h, bcs=bc,
      solver_parameters={
          'ksp_type': 'cg',
          'pc_type': 'gamg',
      },
      options_prefix='poisson')

VTKFile('pvd/poisson_example.pvd').write(u_h)
