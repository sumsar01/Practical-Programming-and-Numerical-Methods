CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

.PHONEY : all
all : linear.out.txt airy.svg
	cat $(word 2,$^)

linear.out.txt :linear
	./$< > $@

airy.svg : airy.gpi airy.data
	gnuplot $<

airy.data : airy
	./$< > $@



.PHONEY : clean
clean :
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
