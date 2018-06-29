CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
objects = main QR_decomp *.dat *.o *.svg *.txt

.PHONEY: all clean

all: output.txt Makefile error_mc.svg
	cat output.txt
	cat errors.txt

error_mc.svg: error_mc.gpi error_mc.txt
	gnuplot $<

output.txt error_mc.txt: main
	./$< 1> output.txt 2> errors.txt

main: main.o MC_int.o lin_alg_fun.o MC_int.h

lin_alg_fun.o:
main.o:
MC_int.o:

clean:
	$(RM) $(objects)
