# Build script for spe2nc. Reads in spectral initial data and writes a netcdf file
# This Makefile ASSUMES $NETCDF points to where netcdf library is installed.
# There should be enough metadata in the netcdf file to use ferret, ncview, etc.

FC         = ifort
OMPFLAG    = -openmp
FFLAGS     = -g -O3 -I$(NETCDF)/include $(OMPFLAG)
FIXEDFLAGS = -g -O2 -I$(NETCDF)/include
LDFLAGS    = -L$(NETCDF)/lib -lnetcdf $(OMPFLAG)
OBJS       = nc.o sigio_module.o spe2nc.o sptsnv.o sptranf0.o sptranf.o sptranf1.o sptranfv.o spuv2dz.o \
             spdz2uv.o spffte.o spwget.o ncpus.o sptez.o sptezm.o sptezmv.o splat.o spsynth.o spanaly.o \
             splegend.o speps.o

%.o: %.F90
	$(FC) -c $(FFLAGS) $<

%.o: %.f
	$(FC) -c $(FIXEDFLAGS) $<

all: spe2nc speinfo

spe2nc: $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

speinfo: speinfo.o sigio_module.o
	$(FC) -o $@ $? $(OMPFLAG)

clean:
	$(RM) $(OBJS) *.mod spe2nc

nc.o:
sigio_module.o: sigio_module.F90
	$(FC) -c $(FIXEDFLAGS) $<
spe2nc.o: nc.o sptsnv.o
sptranf0.o:
sptranf.o:
sptranf1.o:
sptranfv.o:
spuv2dz.o:
spdz2uv.o:
spffte.o:
spwget.o:
ncpus.o:
sptez.o:
sptezm.o:
sptezmv.o:
splat.o:
spsynth.o:
spanaly.o:
splegend.o:
speps.o:
