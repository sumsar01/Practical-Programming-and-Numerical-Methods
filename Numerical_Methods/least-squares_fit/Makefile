CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
object = main QR_decomp fit *.txt *.o *.svg

.PHONEY: all clean

all: output.txt Makefile fit1.svg fit2.svg
	cat output.txt

fit1.svg: fit_1.gpi fit1.txt
	gnuplot $<

fit2.svg: fit_2.gpi fit2.txt
	gnuplot $<

output.txt fit1.txt fit2.txt: main
	./$< 1> output.txt 2> errors.txt

main: main.o QR_decomp.o fit.o fit.h

main.o:
QR_decomp.o:
fit.o:


clean:
	$(RM) $(object)
