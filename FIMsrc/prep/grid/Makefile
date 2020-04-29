# grid Makefile

include ../../macros.make

CNTLDIR      = ../../cntl
CNTLOBJS     = $(CNTLDIR)/fimnamelist.o
GETLVL       = $(BINDIR)/getlvl
GINFO        = $(BINDIR)/ginfo
GINFO_OBJS   = datastru.o GetRegions.o GridGen.o hilbert.o ijblock.o IJblockLayout.o mpi_stubs.o perm.o SquareDecomp.o SquareLayout.o 
GRID         = $(BINDIR)/grid
GRID_OBJS    = avg.o bisect_triangle.o datastru.o icos.o middle.o rotate.o third.o top_gridpoints.o top_triangles.o trisect.o trisect_triangle.o 
INCS         = -I../../cntl -I$(UTILDIR) -I../incmod -I../lsfit
LIBCNTL      = $(LIBDIR)/libcntl.a
LIBSYSSHARE  = $(LIBDIR)/libsysshare.a
LL2IPNPROG   = ll2ipn
SLINTLIB     = $(LIBDIR)/libslint.a
LIBLSFIT     = ../lsfit/liblsfit.a
UTILDIR      = ../../utils
UTILOBJS     = $(UTILDIR)/headers.o $(UTILDIR)/module_initial_chem_namelists.o
#MKLIB = -lmkl_intel_lp64 -lmkl_core -lmkl_sequential

override FFLAGS += $(INCS) $(MICFLAGS)

ifeq ($(P),p)
  SMSLIB = -L$(SMS)/lib -lsms
else
  CPPFLAGS += -DSERIAL
endif

.SUFFIXES:
.PRECIOUS: %.f90

%.f90: %.F90
	$(CPP) $(CPPFLAGS) $< > $@

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

all: $(GRID) $(GINFO) $(GETLVL) $(LL2IPNPROG)

$(GRID) $(GINFO) $(GETLVL): DEPENDENCIES

$(GRID): $(LIBCNTL)(module_control.o) $(LIBSYSSHARE)(sys_share.o) $(GRID_OBJS) $(UTILOBJS) $(CNTLOBJS)
	$(FC) $(FFLAGS) -o $(GRID) $(GRID_OBJS) $(UTILOBJS) $(CNTLOBJS) $(LIBCNTL) $(LIBSYSSHARE) $(LDFLAGS)

$(GINFO): $(LIBCNTL)(module_control.o) $(LIBSYSSHARE)(sys_share.o) $(GINFO_OBJS) $(UTILOBJS) $(CNTLOBJS)
	$(FC) $(FFLAGS) -o $(GINFO) $(GINFO_OBJS) $(LIBCNTL) $(UTILOBJS) $(SMSLIB) $(LDFLAGS) $(LIBSYSSHARE)

$(GETLVL): getlvl.o  $(LIBCNTL)(module_control.o) $(LIBSYSSHARE)(sys_share.o) $(UTILOBJS) $(LIBLSFIT) 
	$(FC) $(FFLAGS) -o $(GETLVL) $< $(LIBCNTL) $(LIBSYSSHARE) $(UTILOBJS) $(LIBLSFIT) $(MKLIB) 

$(LL2IPNPROG): ll2ipn.o
	$(FC) $(FFLAGS) -o $(LL2IPNPROG) $< $(UTILOBJS) $(SLINTLIB)

DEPENDENCIES:
	$(RM) -f Filepath Srcfiles
	echo "." > Filepath
	ls *.F90 > Srcfiles
	$(MKDEPENDS) -m Filepath Srcfiles > $@

-include DEPENDENCIES

clean:
	$(RM) *.o *.mod DEPENDENCIES Filepath Srcfiles
