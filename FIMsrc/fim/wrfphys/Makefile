# wrfphys Makefile

include ../../macros.make
include FIM_WRFP_OBJECTS

INCS = -I../../cntl -I../column -I../../utils

%.o: %.F
	$(FC) -c $(FFLAGS) $(INCS) $(FREEFLAG) $(OPTFLAGS) $<

all: DEPENDENCIES $(OBJS)

-include DEPENDENCIES

DEPENDENCIES: $(wildcard *.F)
	$(RM) Filepath Srcfiles
	echo "." > Filepath
	ls *.F > Srcfiles
	$(MKDEPENDS) -m Filepath Srcfiles > $@

clean:
	$(RM) *.o *.mod *.a DEPENDENCIES Filepath Srcfiles
