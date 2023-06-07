# Jupyter-book


## Pygments

### Lexers
https://pygments.org/docs/lexers/


### Extension for pygments

文件列表

```console
$ tree -L 2 --filesfirst --charset=ascii ../extension 
../extension
|-- README.md
|-- pyproject.toml
`-- xyzutils
    |-- __init__.py
    |-- colab_python.py
    `-- __pycache__

2 directories, 4 files
```

:::{literalinclude} ../extension/pyproject.toml
:caption:
:::

:::{literalinclude} ../extension/xyzutils/__init__.py
:caption:
:::

:::{literalinclude} ../extension/xyzutils/colab_python.py
:caption:
:::