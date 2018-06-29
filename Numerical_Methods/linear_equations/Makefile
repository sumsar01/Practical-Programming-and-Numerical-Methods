CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all clean

all: 		output.txt
output.txt:	main;				./main 1> output.txt 2> errors.txt
main:		main.o test.o lin_eq.o
main.o:
lin_eq.o:
test.o:

clean:
	$(RM) main *.txt *.o
