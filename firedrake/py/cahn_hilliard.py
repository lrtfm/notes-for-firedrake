# split method for Cahn-Hilliard equation
#
# Date: 2023-05-06

from firedrake import *
from firedrake.petsc import PETSc

class Bar(ProgressBar):
    suffix = '%(index)s/%(max)s [%(elapsed_td)s/%(eta_td)s]'
    bar_prefix = ' |'
    bar_suffix = '| '
    empty_fill = ' '
    fill = '#'
    color = None

def u0(x, y):
    return 0.05*cos(2*pi*x)*cos(2*pi*y)

def f_plus(u):
    return u**3

def f_minus(u):
    return u

opts = PETSc.Options()
degree = opts.getInt('degree', default=1)
N = opts.getInt('N', default=100)
M = opts.getInt('M', default=1600)
tau = opts.getReal('tau', default=1e-4)
epsilon = opts.getReal('epsilon', default=0.05)
periodic = opts.getBool('periodic', default=True)

dt = Constant(tau)

if periodic:
    filename = 'pvd/test_ch_periodic.pvd'
    mesh = PeriodicRectangleMesh(N, N, 2, 2)
else:
    filename = 'pvd/test_ch_neumann.pvd'
    mesh = RectangleMesh(N, N, 2, 2)

mesh.coordinates.dat.data[:] = mesh.coordinates.dat.data_ro - 1

V = FunctionSpace(mesh, 'CG', degree)
W = V*V
v, v_test = Function(W), TestFunction(W)
u, w = split(v)
u_test, w_test = split(v_test)

vn = Function(W)
un, wn = vn.subfunctions
un.rename('u')
wn.rename('w')

a = 1/dt*inner(u - un, u_test)*dx + inner(grad(w), grad(u_test))*dx \
    + inner(w, w_test)*dx - epsilon**2*inner(grad(u), grad(w_test))*dx \
    - inner(f_plus(u) - f_minus(un), w_test)*dx

prob = NonlinearVariationalProblem(a, v)
solver = NonlinearVariationalSolver(prob, options_prefix="ch") #  solver_parameters={'snes_monitor': None, 'snes_view': None})

t = 0
x, y = SpatialCoordinate(mesh)
un.interpolate(u0(x,y))

PETSc.Sys.Print(f'Will save result in {filename}')
output = File(filename)
output.write(un, wn, time=t)

for i in Bar('Timestep').iter(range(M)):
    t = (i+1)*tau
    solver.solve()

    vn.assign(v)
    if (i+1)%100 == 0:
        output.write(un, wn, time=t)
