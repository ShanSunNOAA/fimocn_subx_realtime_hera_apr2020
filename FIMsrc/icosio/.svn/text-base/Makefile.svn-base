# icosio Makefile

# CPPFLAGS may include:
#   -DSERIAL for no-MPI build
#   -DNOGRIB for no-GRIB build
#   -DUSEMPIMOD => "use mpi"
#   -UUSEMPIMOD => "include 'mpif.h'"

.PRECIOUS: %.f90

%.f90: %.F90
	$(CPP) $(CPPFLAGS) $< > $@

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

all: icosio.o

icosio.f90: icosio_in_receive.decl icosio_in_receive.exec
icosio.o: ignore_unused.o stringlist.o

clean:
	$(RM) *.f90 *.o *.mod
