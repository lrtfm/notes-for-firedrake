.PHONY : all html pdf clean

all: local

OUTPUT=${HOME}/.aux/firedrake-notes

pdf:
	jupyter-book build --path-output ${OUTPUT}  --builder pdflatex ./
	cp ${OUTPUT}/_build/latex/firedrake-notes.pdf ./

html:
	jupyter-book build --path-output ${OUTPUT} ./

push: html
	@echo Syncing to zzyang.net
	@rsync -rP --delete ${OUTPUT}/_build/html/ zzyang.net:/var/www/html/firedrake-notes/ > rsync.logs
	@echo You can look at your book by click the link
	@echo "\n\thttp://zzyang.net/firedrake-notes/index.html"

local: html
	@echo Syncing to local directory
	@rsync -rP --delete ${OUTPUT}/_build/html/ ~/opt/Sites/firedrake-notes/ > rsync-local.logs
	@echo "\nYou can look at your book by click the link"
	@echo "\n\thttp://localhost/~zzyang/firedrake-notes/index.html"

clean:
	jupyter-book clean ./
	rm -rf ${OUTPUT}/_build/html/
	rm -rf ${OUTPUT}/_build/latex/

