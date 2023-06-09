.PHONY : all html pdf clean

all: html

pdf:
	jupyter-book build --builder pdflatex ./
	cp _build/latex/firedrake-notes.pdf ./

html:
	jupyter-book build ./

clean:
	jupyter-book clean ./

