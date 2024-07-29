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

<!--
之前我一直使用 MATLAB 做有限元计算. 由于它是一个通用编程脚本, 实现复杂有限元程序需要做大量工作,
另外也是商业软件, 我希望找到一个易用的有限元工具包.
最先找到的是 PETSc, 只不过它内容有限元框架开发时间短, 用户较少. 后来就看到了基于 PETSc 的 Firedrake.
Firedrake 是一个求解偏微分方程的有限元框架, 和 FEniCS 有很多类似之处. 不过它主要是由 Python 编写的.
对它内部做一些适应我们算法的修改相对容易, 因此我决定开始学习和使用 Firedrake.
-->

由于 Firedrake 处于开发中, 尚未发布正式版, 接口也时有变动, 给学习和使用带来了一些困难.
在这种情况下, 主要的学习方法是结合它的文档和代码一起看, 有问题时在社区提问.

在学习和使用中, 我遇到过很多问题. 为了便于后续查找, 写了一些笔记, 目前很多内容过于简略.
不过第一二章稍微详细一些, 可以作为初学者的参考. 计划最近对笔记内容进行更新和补充, 使其更加完整.