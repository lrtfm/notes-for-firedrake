
BASE_NAME=01_firedrake_install 02_firedrake_notes

.PHONY : all pdf tex
all: pdf
pdf: $(addsuffix .pdf, $(BASE_NAME))
tex: $(addsuffix .tex, $(BASE_NAME))

template_arg=--TemplateExporter.extra_template_basedirs=./templates
%.pdf: %.ipynb
	jupyter nbconvert $(template_arg) --to pdf $<

%.tex: %.ipynb
	jupyter nbconvert $(template_arg) --to latex $<

.PHONY : test
test:
	echo $(addsuffix .pdf, $(BASE_NAME))

.PHONY : clean
clean:
	$(foreach name,$(BASE_NAME),rm -rf $(name)_files $(name).tex $(name).pdf;)
