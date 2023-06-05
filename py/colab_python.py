from pygments.lexer import RegexLexer, inherit
from pygments.lexers import PythonLexer
from pygments.token import *


class ColabPythonLexer(PythonLexer):
    """All your lexer code goes here!"""
    name = 'colabPython'
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