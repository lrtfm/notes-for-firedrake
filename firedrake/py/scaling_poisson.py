"""用 PETSc 自定义 Stage 划分 Poisson 求解的各步, 用于强扩展性测试.

在 firedrake 目录下运行, 例如:

    mpiexec -n 4 python3 py/scaling_poisson.py -log_view :poisson-np4.csv:ascii_csv

可用的命令行参数:

    -N <int>       粗网格每个方向的单元数 (默认 128)
    -refine <int>  在粗网格上加密的次数, 同时构成多重网格的层次 (默认 3)

各 Stage 的耗时在 -log_view 输出的 "Summary of Stages" 中, 用 ascii_csv 格式时
每个 Stage 有一行事件名为 summary 的记录.
"""
from firedrake import *
from firedrake.petsc import PETSc

opts = PETSc.Options()
N = opts.getInt('N', 128)
refine = opts.getInt('refine', 3)

with PETSc.Log.Stage('MyMesh'):
    coarse = UnitSquareMesh(N, N)
    mesh = MeshHierarchy(coarse, refine)[-1]

with PETSc.Log.Stage('MySpace'):
    V = FunctionSpace(mesh, 'CG', 1)
    x, y = SpatialCoordinate(mesh)
    u, v = TrialFunction(V), TestFunction(V)
    a = inner(grad(u), grad(v))*dx
    L = inner(sin(pi*x)*sin(pi*y), v)*dx
    bc = DirichletBC(V, 0, 'on_boundary')
    u_h = Function(V, name='u_h')

params = {
    'ksp_type': 'cg',
    'ksp_rtol': 1e-8,
    'pc_type': 'mg',              # 几何多重网格, 层次来自 MeshHierarchy
}
problem = LinearVariationalProblem(a, L, u_h, bcs=bc)
solver = LinearVariationalSolver(problem, solver_parameters=params)

with PETSc.Log.Stage('MySolve'):
    solver.solve()

its = solver.snes.ksp.getIterationNumber()
PETSc.Sys.Print(f'nproc = {mesh.comm.size}, dofs = {V.dim()}, CG iterations = {its}')
