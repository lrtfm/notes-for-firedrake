"""Drop English-only content from the LaTeX/PDF build.

The HTML site offers a zh/en toggle (see _static/lang-toggle.js): nodes whose
classes contain ``lang-en`` are shown only in English mode. LaTeX cannot
toggle, so this extension removes those nodes from the doctree when building
with the latex builder, keeping the PDF Chinese-only. Nodes marked
``lang-zh`` are kept (their class is meaningless in LaTeX and ignored).
"""


def _has_lang_en(node):
    try:
        classes = node.get("classes") or []
    except AttributeError:
        return False
    return "lang-en" in classes


def on_doctree_resolved(app, doctree, docname):
    if app.builder.format != "latex":
        return
    for node in list(doctree.findall(_has_lang_en)):
        if node.parent is not None:
            node.parent.remove(node)


def setup(app):
    app.connect("doctree-resolved", on_doctree_resolved)
    return {
        "version": "0.1",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
