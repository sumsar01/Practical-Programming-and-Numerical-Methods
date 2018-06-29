CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
object = *.o *.txt *.svg main *.pdf main jacobi_sweep_timer jacobi_diag_ev_by_ev_one_time jacobi_diag_ev_by_ev_full_time



.PHONEY: all clean

all: output.txt Makefile time.svg
	cat output.txt

output.txt: main
	./$< 1> output.txt

main: main.o print.o jacobi.h jacobi_diag.o

time.svg: time_sweep.txt time_ev_by_ev_one.txt time_ev_by_ev_full.txt
	gnuplot time.gpi

time_sweep.txt: jacobi_sweep_timer
	n=100; $(RM) $@; \
	while [ $${n} -le 600 ] ; do \
		\time --format "%U" --append --output $@ ./$< $$n > /dev/null; \
		n=`expr $$n + 30`; \
	done;

time_ev_by_ev_one.txt: jacobi_diag_ev_by_ev_one_time
	n=100; $(RM) $@; \
	while [ $${n} -le 600 ] ; do \
		\time --format "%U" --append --output $@ ./$< $$n > /dev/null; \
		n=`expr $$n + 30`; \
	done;

time_ev_by_ev_full.txt: jacobi_diag_ev_by_ev_full_time
	n=100; $(RM) $@; \
	while [ $${n} -le 600 ] ; do \
		\time --format "%U" --append --output $@ ./$< $$n > /dev/null; \
		n=`expr $$n + 30`; \
	done;

jacobi_sweep_timer: jacobi_diag.o jacobi.h print.o jacobi_sweep_timer.o

jacobi_diag_ev_by_ev_full_time: jacobi_diag.o jacobi.h print.o jacobi_diag_ev_by_ev_full_time.o

jacobi_diag_ev_by_ev_one_time: jacobi_diag.o jacobi.h print.o jacobi_diag_ev_by_ev_one_time.o

jacobi_sweep_timer.o:
jacobi_diag_ev_by_ev_full_time.o:
jacobi_diag_ev_by_ev_one_time.o:
main.o:
print.o:
jacobi_diag.o:





clean:
	$(RM) $(object)

