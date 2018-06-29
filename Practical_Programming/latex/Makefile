CDLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)


.PHONEY: all clean

report.pdf: report.tex plot.pdf
	pdflatex $<
	pdflatex $<

plot.pdf: data plot.gpi
	gnuplot plot.gpi

data: main
	./$< -3 3 0.2 > $@

main: main.o

.PHONEY: clean

clean:
	$(RM) main *.dat *.o *.svg *.png *.pdf document.dvi document.ps document.aux document.log
