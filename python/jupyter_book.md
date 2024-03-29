# Jupyter-book

## Install

```console
pip install -U jupyter-book
pip install Pillow
sudo apt install plantuml                 # for sphinxcontrib-plantuml
pip install -U sphinxcontrib-plantuml
pip install -U sphinxcontrib-tikz
pip3 install --upgrade docutils==0.17.1   # https://sourceforge.net/p/docutils/patches/195/
```

## MyST

Directives - a block-level extension point

https://myst-parser.readthedocs.io/en/latest/syntax/roles-and-directives.html#syntax-directives

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