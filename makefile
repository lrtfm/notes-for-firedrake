.PHONY : all html pdf clean

all: html pdf

pdf:
	jupyter-book build --builder pdflatex ./

html:
	jupyter-book build ./

clean:
	jupyter-book clean ./

