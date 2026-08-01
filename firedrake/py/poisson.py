"""并行求解 Poisson 方程的最小示例.

在 firedrake 目录下运行:

    mpiexec -n 2 python3 py/poisson.py
"""
from firedrake import *
from firedrake.petsc import PETSc

N = 4
test_mesh = RectangleMesh(N, N, 1, 1)
x, y = SpatialCoordinate(test_mesh)
f = sin(pi*x)*sin(pi*y)
g = Constant(0)

V = FunctionSpace(test_mesh, 'CG', degree=1)

u, v = TrialFunction(V), TestFunction(V)
a = inner(grad(u), grad(v))*dx
L = inner(f, v)*dx

bc = DirichletBC(V, g=g, sub_domain='on_boundary')

u_h = Function(V, name='u_h')
solve(a == L, u_h, bcs=bc)

# assemble 是集合操作, 所有进程一起调用; Print 只由 0 号进程输出
PETSc.Sys.Print(f'[{COMM_WORLD.size} process(es)] integral of u_h: '
                f'{assemble(u_h*dx):.10e}')

output = VTKFile('data/result.pvd')
output.write(u_h)
