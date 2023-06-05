# Try Firedrake on Colab

Colab is short for Colaboratory (which can be considered as an online version of Jupyter, allowing you to write and execute Python code in your browser).

[FEM on Colab](https://fem-on-colab.github.io/index.html) supports the installation of [FEniCS, FEniCSx, Firedrake, NGSolve, gmsh](https://fem-on-colab.github.io/packages.html) on colab
## Import package

### Firedrake

About 3 minutes

:::python
try:
    import firedrake
except ImportError:
    !wget "https://fem-on-colab.github.io/releases/firedrake-install-real.sh" \
        -O "/tmp/firedrake-install.sh" && bash "/tmp/firedrake-install.sh"
    import firedrake
:::

### Gmsh

:::python
try:
    import gmsh
except ImportError:
    !wget "https://fem-on-colab.github.io/releases/gmsh-install.sh" \
        -O "/tmp/gmsh-install.sh" && bash "/tmp/gmsh-install.sh"
    import gmsh
:::

## Examples

https://colab.research.google.com/drive/1gM3zMWTskH7XyDi1yJL76BPFnOJjSdYh?usp=sharing
