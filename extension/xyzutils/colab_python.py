from pygments.lexer import RegexLexer, inherit
from pygments.lexers import PythonLexer
from pygments.token import *

__all__ = ( "ColabPythonLexer", )

class ColabPythonLexer(PythonLexer):
    """All your lexer code goes here!"""
    name = 'colab-python'
    aliases = ['colabpy']
    filenames = ['*.colabpy']

    tokens = {
        'root': [
            (r'\s+\!([^\n]+\\\n)*[^\n]*[^\\]\n', String),
            inherit
        ]
    }

if __name__ == '__main__':
    pass 