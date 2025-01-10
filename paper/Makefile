SRC := paper.tex \
       abstract.tex \
	   intro.tex \
       preliminaries.tex \
	   accuracy.tex \
	   subset.tex \
	   related.tex \
	   evaluation.tex \
	   conclusion.tex \
       paper.bib

all: paper.pdf outline.pdf

paper.pdf: $(SRC)
	latexmk -pdf $<

outline.pdf: outline.tex
	latexmk -pdf $<

paper.zip: $(SRC) paper.bbl lipics-v2021.cls
	$(RM) paper.zip
	zip paper.zip $(SRC) latex/model.tex paper.bbl lipics-v2021.cls

clean:
	$(RM) *.aux *.log *.out *.vrb paper.pdf \
              *.bbl *.blg *.fdb_latexmk *.toc *.fls *.cut
