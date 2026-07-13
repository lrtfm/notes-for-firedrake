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

:::lang-zh
[Firedrake](https://www.firedrakeproject.org/) 是一个用于求解偏微分方程的有限元框架,
类似于 FEniCS, 且主要由 Python 编写, 方便根据算法需求进行修改.

在使用 Firedrake 的过程中, 我记录了一些学习笔记, 以便后续参考.
当前笔记内容较为简略, 但第一二章较为详细, 可供初学者参考.
本笔记以中文为主, 正在逐步补充英文内容:
点击页面右上角的 **EN** 按钮可切换到英文模式.
:::

:::lang-en
[Firedrake](https://www.firedrakeproject.org/) is a finite element framework
for solving partial differential equations, similar to FEniCS.
It is written mainly in Python, which makes it easy to adapt to the needs of
your own algorithms.

These are notes I took while learning and using Firedrake, kept here for
future reference. Most of the notes are brief, but the first two chapters are
fairly detailed and should be helpful for beginners. The notes are written
mostly in Chinese and English content is being added gradually: untranslated
parts are shown in Chinese. Click the **中文** button at the top of the page
to switch back to Chinese.
:::

## [阅读指南]{.lang-zh} [How to read]{.lang-en}

:::lang-zh
+ 正文分为三个部分: **入门 (Tutorials)** 介绍求解各类方程的基本流程;
  **进阶专题 (Advanced Topics)** 包括并行计算, PETSc, 网格生成, 调试与性能分析;
  **内部机制 (Internals)** 探讨 Firedrake 的实现细节.
  初学者建议从 [Poisson 方程 I](firedrake/poisson_I.ipynb) 和
  [Poisson 方程 II](firedrake/poisson_II.ipynb) 开始.
  每章开头附有英文概要 (Overview), 并注明运行前提.

+ 安装方法请参考 Installation 部分.
  推荐初学者使用 [Docker](https://www.firedrakeproject.org/install.html#docker) 或
  [Colab](installation/colab.md), 无需本地编译.

+ 附加信息部分收录了 Python, Linux, Git 和编辑器等相关工具的笔记,
  供不熟悉这些工具的读者参考.

+ 常用术语的中英文对照见[术语对照表](reference/glossary.md).
:::

:::lang-en
+ The main text has three parts: *Tutorials* on the basic workflow of solving
  PDEs, *Advanced Topics* (parallel computing, PETSc, meshing, debugging and
  profiling), and *Internals* on implementation details. Beginners are
  encouraged to start with the two chapters on the Poisson equation. Each
  chapter starts with a short English overview stating its scope and
  prerequisites.

+ See the Installation part for how to install Firedrake. Beginners may find
  [Docker](https://www.firedrakeproject.org/install.html#docker) or
  [Colab](installation/colab.md) the easiest way to get started, as no local
  build is needed.

+ The Additional Information part collects notes on related tools such as
  Python, Linux, Git and editors, for readers who are not yet familiar with
  them.

+ See the [Glossary](reference/glossary.md) for Chinese–English terminology
  used throughout the notes.
:::

## [本地构建]{.lang-zh} [Build locally]{.lang-en}

:::lang-zh
本书使用 [Jupyter Book](https://jupyterbook.org/) 构建:
:::

:::lang-en
The book is built with [Jupyter Book](https://jupyterbook.org/):
:::

```bash
python3 -m pip install -r requirements.txt
jupyter-book build ./
```

:::lang-zh
CI 使用最新的 `firedrakeproject/firedrake-vanilla-default` 镜像构建.
推送前可运行 `make docker-html`, 在与 CI 相同的镜像中做一次验证构建
(先用 `make docker-pull` 更新镜像), 以提前发现 Firedrake 新版本的接口变化.
:::

:::lang-en
The CI builds with the latest `firedrakeproject/firedrake-vanilla-default`
image. Before pushing, you can run `make docker-html` to verify the build
inside the same image (update it first with `make docker-pull`), catching
API changes of new Firedrake releases early.
:::

```{nb-exec-table}
```
