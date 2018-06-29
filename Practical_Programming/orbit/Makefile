CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
PLOTTER = gnuplot

PLOTS = orbit1.svg orbit2.svg orbit3.svg

.PHONEY: default
default: $(PLOTS) log.svg

log.svg: log.gpi log.dat
	gnuplot $<

log.dat: log_function
	./$< > $@

 $(PLOTS): orbit.gpi orbit.dat
	gnuplot $<

orbit.dat: main
	./$< > $@

#out.txt : main orbit.gpi 
#	; ./$< > $@ && gnuplot orbit.gpi


.PHONEY: clean
clean: ; $(RM) $(PLOTS) main log.dat log.svg


