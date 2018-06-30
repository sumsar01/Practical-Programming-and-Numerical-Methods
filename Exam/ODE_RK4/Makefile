CDLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
objects = main QR_decomp *.dat *.o *.svg *.txt


.PHONEY: all clean


all: output.txt Makefile airy.svg sin.svg
	cat output.txt
	cat errors.txt

airy.svg: airy.gpi airy.txt
	gnuplot airy.gpi

sin.svg: sin.gpi sin.txt
	gnuplot sin.gpi

output.txt airy.txt sin.txt: main
	./$< 1> output.txt 2> errors.txt

main: main.o QR_decomp.o ODE.o lin_alg_fun.o ODE.h

lin_alg_fun.o:
main.o:
QR_decomp.o:
ODE.o:

clean:
	$(RM) $(objects)
