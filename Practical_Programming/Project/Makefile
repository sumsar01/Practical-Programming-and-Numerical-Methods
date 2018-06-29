CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all clean

all: project.tex plot.pdf
	pdflatex $<
	pdflatex $<

plot.pdf: plot.gpi data
	gnuplot $<

data: main
	./$< > $@

main: main.o root.o

clean:
	$(RM) main *.o *pdf data document.dvi document.ps document.aux document.log
