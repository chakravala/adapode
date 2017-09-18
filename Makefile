# adapode, Copyright (C) 2013 Michael Reed

FCOMPLR=gfortran -g
CFLAGS=-O1
OBJECTS=Constants.o Functions.o OdeEuler.o OdeRungeKutta.o OdeMultistep.o OdeSolve.o OdePrograms.o adapode.o

adapode: $(OBJECTS)
	$(FCOMPLR) $(OBJECTS) -o $@
	- rm -f *.o *~ *.mod

%.o: %.f95
	$(FCOMPLR) $(CFLAGS) -c $<

run:
	./adapode
	gnuplot> load 'Lorenz.plt'
	- rm -f load

clean:
	- rm -f *.o *~ *.mod

cleanall:
	- rm -f *.out fort.* *.err core* *.o *~ *.mod adapode d load
