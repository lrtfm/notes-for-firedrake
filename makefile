.PHONY : all html pdf clean

all: html

OUTPUT=${HOME}/aux/firedrake

pdf:
	jupyter-book build --path-output ${OUTPUT}  --builder pdflatex ./
	cp ${OUTPUT}/_build/latex/firedrake-notes.pdf ./

html:
	jupyter-book build --path-output ${OUTPUT} ./

push: html
	rsync -rP --delete ${OUTPUT}/_build/html/ zzyang.net:/var/www/html/firedrake-notes/ > rsync.logs

clean:
	jupyter-book clean ./
	rm -rf ${OUTPUT}/_build/html/
	rm -rf ${OUTPUT}/_build/latex/

