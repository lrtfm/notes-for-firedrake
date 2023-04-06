from firedrake import *
from firedrake.petsc import PETSc

rank, size = COMM_WORLD.rank, COMM_WORLD.size

opts = PETSc.Options()
N = opts.getInt('N', 32*size)

test_mesh = RectangleMesh(nx=N, ny=N, Lx=1, Ly=1)   # Build mesh on the domain
x, y = SpatialCoordinate(test_mesh)
f = sin(pi*x)*sin(pi*y)
g = Constant(0)

V = FunctionSpace(test_mesh, 'CG', degree=1)        # define FE space

u, v = TrialFunction(V), TestFunction(V)            # define trial and test function 
a = inner(grad(u), grad(v))*dx
L = inner(f, v)*dx                                  # or f*v*dx

bc = DirichletBC(V, g=g, sub_domain='on_boundary')

u_h = Function(V, name='u_h')
solve(a == L, u_h, bcs=bc, options_prefix='')          # We will introduce other ways to code in the following.

