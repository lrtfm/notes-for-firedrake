---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Firedrake Notes

```{raw} latex
\addcontentsline{toc}{chapter}{Preface}
\chapter*{Preface}
\markboth{Preface}{Preface}
```

我曾主要使用 MATLAB 进行有限元计算.
由于 MATLAB 是一种通用编程语言, 实现复杂的有限元程序需要投入大量额外工作.
另外其为商业软件, 使用上有一些限制, 因此我希望找到一个更易用且开源的有限元工具包.

最初接触的是 PETSc, 然而由于其有限元框架开发时间较短, 用户较少, 相关资源也较为有限.
随后, 我发现了基于 PETSc 的 Firedrake. Firedrake 是一个用于求解偏微分方程的有限元框架,
类似于 FEniCS, 且主要由 Python 编写, 方便根据算法需求进行修改.

由于 Firedrake 仍在开发中, 学习方式主要依赖于文档与源码的结合, 并通过社区寻求解决方案.

在使用 Firedrake 的过程中, 我记录了一些学习笔记, 以便后续参考.
当前笔记内容较为简略, 但第一二章较为详细, 可供初学者参考.

<!-- 我计划近期对笔记进行更新和完善，使其内容更为系统和全面. -->