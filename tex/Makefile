all:report.pdf #preview

report.pdf:report.tex ./sec/*.tex
	xelatex report.tex
	xelatex report.tex

preview:report.pdf
	evince report.pdf&

.PHONY:clean realclean
clean:
	-rm -f *.aux *.dvi *.log *.out *.blg *.bbl *.out

realclean:
	-rm -f report.pdf
