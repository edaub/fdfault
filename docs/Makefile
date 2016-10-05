## userguide.pdf : Create PDF documentation (same as make with no arguments)
userguide.pdf: userguide.tex
	pdflatex userguide.tex
	pdflatex userguide.tex
	rm -f userguide.aux userguide.log userguide.toc

.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<