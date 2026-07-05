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

[Firedrake](https://www.firedrakeproject.org/) 是一个用于求解偏微分方程的有限元框架,
类似于 FEniCS, 且主要由 Python 编写, 方便根据算法需求进行修改.

在使用 Firedrake 的过程中, 我记录了一些学习笔记, 以便后续参考.
当前笔记内容较为简略, 但第一二章较为详细, 可供初学者参考.
本笔记目前以中文为主, 部分章节为英文, 后续会逐步补充双语内容.

[Firedrake](https://www.firedrakeproject.org/) is a finite element framework
for solving partial differential equations, similar to FEniCS.
It is written mainly in Python, which makes it easy to adapt to the needs of
your own algorithms.

These are notes I took while learning and using Firedrake, kept here for
future reference. Most of the notes are brief, but the first two chapters are
fairly detailed and should be helpful for beginners. The notes are currently
written mostly in Chinese, with some sections in English; bilingual content
will be added gradually.

## 阅读指南 (How to read)

+ 正文分为三个部分: **入门 (Tutorials)** 介绍求解各类方程的基本流程;
  **进阶专题 (Advanced Topics)** 包括并行计算, PETSc, 网格生成, 调试与性能分析;
  **内部机制 (Internals)** 探讨 Firedrake 的实现细节.
  初学者建议从 [Poisson 方程 I](firedrake/poisson_I.ipynb) 和
  [Poisson 方程 II](firedrake/poisson_II.ipynb) 开始.
  每章开头附有英文概要 (Overview), 并注明运行前提.

  The main text has three parts: *Tutorials* on the basic workflow of solving
  PDEs, *Advanced Topics* (parallel computing, PETSc, meshing, debugging and
  profiling), and *Internals* on implementation details. Beginners are
  encouraged to start with the two chapters on the Poisson equation. Each
  chapter starts with a short English overview.

+ 常用术语的中英文对照见[术语对照表](reference/glossary.md).

  See the [Glossary](reference/glossary.md) for Chinese–English terminology.

+ 安装方法请参考 [Installation](installation/install.md) 部分.
  推荐初学者使用 [Docker](installation/docker.ipynb) 或
  [Colab](installation/colab.md), 无需本地编译.

  See the Installation part for how to install Firedrake. Beginners may find
  Docker or Colab the easiest way to get started, as no local build is needed.

+ 附加信息部分收录了 Python, Linux, Git 和编辑器等相关工具的笔记,
  供不熟悉这些工具的读者参考.

  The Additional Information part collects notes on related tools such as
  Python, Linux, Git and editors, for readers who are not yet familiar with
  them.

## 本地构建 (Build locally)

本书使用 [Jupyter Book](https://jupyterbook.org/) 构建.
The book is built with Jupyter Book:

```bash
python3 -m pip install -r requirements.txt
jupyter-book build ./
```

```{nb-exec-table}
```
