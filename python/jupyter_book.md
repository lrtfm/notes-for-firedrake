# Jupyter-book

## Install

```console
sudo apt install plantuml                 # for sphinxcontrib-plantuml
```

```
pip install -U jupyter-book
pip install --upgrade docutils
pip install Pillow
pip install -U sphinxcontrib-plantuml
pip install -U sphinxcontrib-tikz
```

## MyST

Directives - a block-level extension point

https://myst-parser.readthedocs.io/en/latest/syntax/roles-and-directives.html#syntax-directives

## 双语内容写法 (Bilingual content)

本书通过页面右上角的按钮切换中英文, 由 [`_static/lang-toggle.js`](../_static/lang-toggle.js)
和 [`_static/lang-toggle.css`](../_static/lang-toggle.css) 实现. 写作约定如下:

+ 未标记的内容 (所有代码单元和未翻译的文字) 在两种语言下都显示;
+ 只在中文模式显示的内容, 用 `lang-zh` 围栏包裹;
+ 只在英文模式显示的内容, 用 `lang-en` 围栏包裹.

在 markdown 文件或 notebook 的 markdown 单元中:

```markdown
:::lang-zh
中文内容.
:::

:::lang-en
English content.
:::
```

若内容中含有其他指令 (嵌套围栏), 外层围栏使用更多冒号 (`::::`).

标题使用行内语法 (需启用 `attrs_inline` 扩展, 已在 `_config.yml` 中开启),
正文标题和侧边栏目录都会随语言切换:

```markdown
# [中文标题]{.lang-zh} [English Title]{.lang-en}
```

对于公式多, 文字少的短句 (一两句话以内), 不必使用成对围栏,
直接用行内语法标记少量文字, 公式只写一遍:

```markdown
[其中]{.lang-zh}[where]{.lang-en} $a > 0$ [且]{.lang-zh}[and]{.lang-en} $b > 0$.
```

注意把空格放在 span 外面, 不要依赖 span 内部的首尾空格.
编号公式 (`\begin{equation}`) 任何时候都不要复制两份, 否则编号会错乱;
将其放在语言块之外共享.

侧边栏的部分标题 (part caption) 来自 `_toc.yml`, 是纯文本,
无法使用上述语法, 由 `lang-toggle.js` 中的 `CAPTIONS` 映射表翻译;
新增或修改 part 时需要同步更新该映射表.
对于 admonition 等指令, 也可以直接添加类名, 例如各章开头的英文概要:

````markdown
```{admonition} 概要 (Overview)
:class: tip lang-en
English overview of this chapter.
```
````

注意: LaTeX/PDF 构建不执行 JavaScript. 本地扩展
[`_ext/lang_filter.py`](../_ext/lang_filter.py) 会在 LaTeX 构建时删除所有
`lang-en` 节点, 因此 PDF 只包含中文内容.

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