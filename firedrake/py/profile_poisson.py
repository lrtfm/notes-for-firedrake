"""用 PETSc 自定义事件标注 Poisson 求解的各个阶段, 配合 -log_view 使用.

在 firedrake 目录下运行:

    python3 py/profile_poisson.py -log_view
    python3 py/profile_poisson.py -log_view :profile.txt:ascii_flamegraph

脚本中没有调用 PETSc.Log.begin(): 命令行传了 -log_view 时 PETSc 已经自己把日志
开好了 (普通格式开的是扁平日志, :...:ascii_xml 和 :...:ascii_flamegraph 开的是
嵌套日志). 这时再调用一次 Log.begin() 会额外启动一个扁平日志处理器; 用后两种
嵌套格式时, 程序退出会报 "Trying to end paused event, not allowed".
"""
from firedrake import *
from firedrake.petsc import PETSc


@PETSc.Log.EventDecorator('MyAssemble')
def build_system(a, L, bc):
    return assemble(a, bcs=bc), assemble(L, bcs=bc)


N = 64
with PETSc.Log.Event('MyMesh'):
    mesh = UnitSquareMesh(N, N)
    V = FunctionSpace(mesh, 'CG', 1)

x, y = SpatialCoordinate(mesh)
u, v = TrialFunction(V), TestFunction(V)
a = inner(grad(u), grad(v))*dx
L = inner(sin(pi*x)*sin(pi*y), v)*dx
bc = DirichletBC(V, 0, 'on_boundary')
u_h = Function(V, name='u_h')

A, b = build_system(a, L, bc)

with PETSc.Log.Event('MySolve'):
    solve(A, u_h, b, solver_parameters={'ksp_type': 'cg',
                                        'pc_type': 'jacobi',
                                        'ksp_rtol': 1e-10})

PETSc.Sys.Print(f'dofs = {V.dim()}, L2 norm of u_h = {norm(u_h):.6e}')
