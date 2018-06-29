CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
object = main QR_decomp print root_finding *.txt *.o *.svg

.PHONEY: all clean

all: output.txt Makefile fit.svg
	cat output.txt
	cat errors.txt

fit.svg: fit.gpi data.txt
	gnuplot fit.gpi

output.txt data.txt: main
	./$< 1> output.txt 2> errors.txt


main: main.o QR_decomp.o minimization.o lin_alg_fun.o minimization.h

main.o:
QR_decomp.o:
minimization.o:
lin_alg_fun.o:

clean:
	$(RM) $(object)

