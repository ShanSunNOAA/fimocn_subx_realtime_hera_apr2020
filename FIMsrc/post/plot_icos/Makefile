# plot_icos Makefile

FC=ifort
FCncarg=ncargf90

BIN = plt_icos
OBJS = pltsubs.o 

%.o: %.f90
	$(FC) $(FLAGS) -c $<

%.o: %.F90
	$(FC) $(FLAGS) -c $<

all: $(BIN)

$(BIN): $(BIN).F90 $(OBJS)
	$(FCncarg) $(FFLAGS) -o plt_icos $^ 

clean:
	$(RM) *.o $(BIN)

tarfile:
	tar cvf $(BIN).tar *.f90 Makefile
