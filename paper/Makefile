SRC := paper.tex \
       abstract.tex \
       preliminaries.tex \
	   intro.tex \
	   quickhull.tex \
	   implementation.tex \
	   accuracy.tex \
	   partition.tex \
	   comparison.tex \
	   conclusion.tex \
       paper.bib

all: paper.pdf

paper.pdf: $(SRC)
	latexmk -pdf $<
	grep Warning paper.log

paper.zip: $(SRC) paper.bbl acmart.cls
	$(RM) paper.zip
	zip paper.zip $(SRC) latex/model.tex paper.bbl acmart.cls

clean:
	$(RM) *.aux *.log *.out *.vrb paper.pdf \
              *.bbl *.blg *.fdb_latexmk *.toc *.fls *.cut
