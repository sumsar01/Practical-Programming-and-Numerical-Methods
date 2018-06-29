CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
object = main QR_decomp print root_finding *.txt *.o *.svg

.PHONEY: all clean

all: output.txt	Makefile
	cat output.txt
	cat errors.txt

output.txt: main
	./$< 1> output.txt 2> errors.txt

main: main.o QR_decomp.o root_finding.o print.o root_finding.h

main.o:
QR_decomp.o:
print.o:
root_finding.o:

clean:
	$(RM) $(object)

