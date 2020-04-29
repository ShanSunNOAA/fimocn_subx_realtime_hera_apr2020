# column_chem Makefile

include ../../macros.make
include FIM_COLUMNC_OBJECTS

INCS = -I../../cntl -I../column -I../../utils

%.o: %.F90
	$(FC) -c $(FFLAGS) $(INCS) $(OPTFLAGS) $<

all: DEPENDENCIES $(OBJS)

-include DEPENDENCIES

DEPENDENCIES:
	$(RM) Filepath Srcfiles
	echo "." > Filepath
	ls *.F90 > Srcfiles
	$(MKDEPENDS) -m Filepath Srcfiles > $@

clean:
	$(RM) *.o *.mod *.a DEPENDENCIES Filepath Srcfiles
