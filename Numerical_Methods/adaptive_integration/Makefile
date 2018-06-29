CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
objects = main QR_decomp *.dat *.o *.svg *.txt

.PHONEY: all clean

all: output.txt Makefile
	cat output.txt
	cat errors.txt

output.txt: main
	./$< 1> output.txt 2> errors.txt

main: main.o adapt_int.o lin_alg_fun.o adapt_int.h

main.o:
adapt_int.o:
lin_alg_fun.o:

clean:
	$(RM) $(objects)


