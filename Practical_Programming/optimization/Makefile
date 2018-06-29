CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
objects = main *.dat *.o *.pdf opg1.txt

.PHONEY: all clean

all: opg1.txt opg2.txt


opg1.txt: minimum
	./$< > $@

opg2.txt: main activity.gpi
	./$< > $@ && gnuplot activity.gpi


clean:
	$(RM) $(objects)
