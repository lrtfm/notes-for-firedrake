.PHONY : all html pdf clean

all: local

NAME=firedrake-notes
OUTPUT=${HOME}/.aux/firedrake-notes
LATEX_BASE=${OUTPUT}/_build/latex
HTML_BASE=${OUTPUT}/_build/html
TEX=${LATEX_BASE}/${NAME}.tex
PDF=${LATEX_BASE}/${NAME}.pdf

# Define a variable to find all .ipynb files excluding .ipynb_checkpoints
NOTEBOOKS := $(shell find . -name "*.ipynb" -not -path "*.ipynb_checkpoints*" -not -path "*publish*")
CHANGED_NOTEBOOKS := $(shell git diff --name-only -- '*.ipynb')
# the path to executed ipynb files
PUBLISH_DIR := ${OUTPUT}/publish

pdf:
	jupyter-book build --path-output ${OUTPUT}  --builder latex ./
	gsed -i.bak -e 's/\([^a-zA-Z0-9\s;,.]\)\s\+\(\\sphinxstylestrong\)/\1\2/g' ${TEX}
	cd ${LATEX_BASE} && latexmk -xelatex ${NAME}.tex
	cp ${PDF} ./

html:
	OMP_NUM_THREADS=1 jupyter-book build --path-output ${OUTPUT} ./

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

# Rule to execute a notebook only if content has changed
$(PUBLISH_DIR)/%.ipynb: %.ipynb
	@mkdir -p $(dir $(PUBLISH_DIR)/$*)
	@# Calculate new checksum
	@md5sum $< | awk '{ print $$1 }' > $(PUBLISH_DIR)/$*.ipynb.checksum.new
	@# Check if checksum has changed
	@if [ ! -f $(PUBLISH_DIR)/$*.ipynb.checksum ] || ! cmp -s $(PUBLISH_DIR)/$*.ipynb.checksum $(PUBLISH_DIR)/$*.ipynb.checksum.new; then \
		python3 -m nbconvert --execute --to notebook --output-dir=$(PUBLISH_DIR) --output $*.ipynb $<; \
		mv $(PUBLISH_DIR)/$*.ipynb.checksum.new $(PUBLISH_DIR)/$*.ipynb.checksum; \
	else \
		rm $(PUBLISH_DIR)/$*.ipynb.checksum.new; \
	fi

# Target to execute all notebooks
execute: $(NOTEBOOKS:%.ipynb=$(PUBLISH_DIR)/%.ipynb)
	@echo "All notebooks have been executed and saved to the $(PUBLISH_DIR) directory."

%.ipynb.clear: %.ipynb
	python3 -m nbconvert \
		--ClearOutputPreprocessor.enabled=True \
		--ClearMetadataPreprocessor.enabled=True \
		--ClearMetadataPreprocessor.preserve_cell_metadata_mask='[("tags")]' \
		--inplace $<

clearoutput: $(CHANGED_NOTEBOOKS:=.clear)
	@echo "All notebook outputs and metadata have been cleared."

clean:
	jupyter-book clean ./
	rm -rf ${HTML_BASE}/
	rm -rf ${LATEX_BASE}/

clean-all:
	rm -rf ${OUTPUT}/_build