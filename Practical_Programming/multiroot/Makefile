CFLAGS = -Wall  -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all

all: problem1.txt problem2.txt
	cat $<

problem1.txt: rosenbrock
	./$< > $@

problem2.txt: main wave.gpi
	./$< > $@ && gnuplot wave.gpi

main: main.o Fe.o

rosenbrock:

clean:
	$(RM) rosenbrock *.o

