# sys_share Makefile

include ../macros.make

FLAGS       = $(FFLAGS) $(MICFLAGS)
LIBSYSSHARE = $(LIBDIR)/libsysshare.a
OBJ         = sys_share.o

.PHONY: all clean

.PRECIOUS: %.f90

%.f90: %.F90
	$(CPP) $(CPPFLAGS) $< > $@

%.o: %.f90
	$(FC) -c $(FLAGS) $<

all: $(LIBSYSSHARE)

$(LIBSYSSHARE): $(LIBSYSSHARE)($(OBJ))
	$(RM) *.o

clean:
	$(RM) *.o *.mod
