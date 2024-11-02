.PHONY : all html pdf clean

all: local

NAME=firedrake-notes
OUTPUT=${HOME}/.aux/firedrake-notes
LATEX_BASE=${OUTPUT}/_build/latex
HTML_BASE=${OUTPUT}/_build/html
TEX=${LATEX_BASE}/${NAME}.tex
PDF=${LATEX_BASE}/${NAME}.pdf

pdf:
	jupyter-book build --path-output ${OUTPUT}  --builder latex ./
	gsed -i.bak -e 's/\([^a-zA-Z0-9\s;,.]\)\s\+\(\\sphinxstylestrong\)/\1\2/g' ${TEX}
	cd ${LATEX_BASE} && latexmk -xelatex ${NAME}.tex
	cp ${PDF} ./

html:
	jupyter-book build --path-output ${OUTPUT} ./

push: html
	@echo Syncing to zzyang.net
	@rsync -rP --delete ${HTML_BASE}/ zzyang.net:/var/www/html/firedrake-notes/ > ${OUTPUT}/rsync.logs
	@echo You can look at your book by click the link
	@echo "\n\thttp://zzyang.net/firedrake-notes/index.html"

local: html
	@echo Syncing to local directory
	@rsync -rP --delete ${HTML_BASE}/ ~/opt/Sites/firedrake-notes/ > ${OUTPUT}/rsync-local.logs
	@echo "\nYou can look at your book by click the link"
	@echo "\n\thttp://localhost/~zzyang/firedrake-notes/index.html"

clearoutput:
	find . -name "*.ipynb" -not -path "*.ipynb_checkpoints*" | \
		xargs python3 -m nbconvert \
		--ClearOutputPreprocessor.enabled=True \
		--ClearMetadataPreprocessor.enabled=True \
		--ClearMetadataPreprocessor.preserve_cell_metadata_mask='[("tags")]' \
		--inplace

execute:
	find . -name "*.ipynb" -not -path "*.ipynb_checkpoints*" | \
		xargs python3 -m nbconvert --execute --inplace

clean:
	jupyter-book clean ./
	rm -rf ${HTML_BASE}/
	rm -rf ${LATEX_BASE}/

