# [术语对照表]{.lang-zh} [Glossary]{.lang-en}

本表给出笔记中常用术语的中英文对照, 以及 Firedrake/UFL 中对应的对象或函数,
方便对照阅读英文文档.

This table lists the Chinese and English terms used in these notes, together
with the corresponding objects or functions in Firedrake/UFL.

| 中文 | English | Firedrake / UFL |
|------|---------|-----------------|
| 网格 | mesh | `Mesh`, `UnitSquareMesh`, `RectangleMesh`, ... |
| 单元 | cell/element | `mesh.ufl_cell()` |
| 面/边/顶点 | facet/edge/vertex | `FacetNormal`, `ds`, `dS` |
| 函数空间 | function space | `FunctionSpace`, `VectorFunctionSpace` |
| 混合函数空间 | mixed function space | `MixedFunctionSpace`, `V * Q` |
| 有限元 | finite element | `FiniteElement` |
| 试探函数 | trial function | `TrialFunction` |
| 检验函数 | test function | `TestFunction` |
| 弱形式/变分形式 | weak/variational form | `a == L` |
| 双线性形式 | bilinear form | `inner(grad(u), grad(v))*dx` |
| 测度 | measure | `dx` (单元), `ds` (外边界), `dS` (内部面) |
| 边界条件 | boundary condition | `DirichletBC` |
| 组装 | assemble | `assemble` |
| 插值 | interpolation | `interpolate` / `Interpolate` |
| 投影 | projection | `project` |
| 自由度 | degree of freedom (DoF) | `Function.dat` |
| 求解器 | solver | `solve`, `LinearVariationalSolver`, `NonlinearVariationalSolver` |
| 线性求解器/Krylov 方法 | Krylov subspace method | PETSc `KSP` (`ksp_type`) |
| 预条件子 | preconditioner | PETSc `PC` (`pc_type`) |
| 非线性求解器 | nonlinear solver | PETSc `SNES` (`snes_type`) |
| 零空间 | nullspace | `VectorSpaceBasis`, `nullspace=` |
| 间断有限元 | discontinuous Galerkin (DG) | `FunctionSpace(mesh, 'DG', k)` |
| 等参元 | isoparametric element | 高阶坐标 (`mesh.coordinates`) |
| 数值积分公式 | quadrature rule | `quadrature_degree` |
| 外法向 | outward normal | `FacetNormal(mesh)` |
| 收敛阶 | convergence order/rate | `errornorm` + 网格序列 |
| 网格自适应 | mesh adaptivity | `adaptMetric`, `adaptLabel` (DMPlex) |
| 并行计算 | parallel computing | `mpiexec`, `mesh.comm` |
