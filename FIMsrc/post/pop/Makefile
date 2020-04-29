# pop Makefile

include ../../macros.make

DEPLIBS     = $(LIBBACIO) $(LIBW3) $(LIBSLINT) $(LIBVLINT) $(LIBGRIBIO) $(LIBCNTL) $(LIBSYSSHARE)
INCS        = -I ../../cntl -I../../prep/incmod -I $(UTILDIR)
LIBBACIO    = $(LIBDIR)/libbacio_4.a
LIBCNTL     = $(LIBDIR)/libcntl.a
LIBGRIBIO   = $(LIBDIR)/libgribio.a
LIBS        = -L$(LIBDIR) -lbacio_4 -lw3_4 -lslint  -lvlint -lgribio -lw3_4 -lcntl -lsysshare -lip
LIBSLINT    = $(LIBDIR)/libslint.a
LIBSYSSHARE = $(LIBDIR)/libsysshare.a
LIBVLINT    = $(LIBDIR)/libvlint.a
LIBW3       = $(LIBDIR)/libw3_4.a
POP         = $(BINDIR)/pop
POPREADINIT = $(BINDIR)/pop_read_init
SLINT_INIT  = $(BINDIR)/slint_init_save
UTILDIR     = ../../utils
UTILOBJS    = $(UTILDIR)/module_initial_chem_namelists.o $(UTILDIR)/headers.o
CNTLDIR     = ../../cntl
CNTLOBJS    = $(CNTLDIR)/module_control.o $(CNTLDIR)/fimnamelist.o

.PRECIOUS: %.f90

.PHONY: all clean

#JR If OPTFLAGS not empty, hopefully compiler will override contradictory FFLAGS settings
%.o: %.F90
	$(FC) -c -g $(FFLAGS) $(OPTFLAGS) $(MICFLAGS) $(INCS) $< 

all: $(POP) $(POPREADINIT) $(SLINT_INIT) post.o 

post.o: post.F90 varinfo.o

$(POP): pop.f90 smooth.o varinfo.o get_tgtgrid.o gaulat.o $(DEPLIBS) $(UTILOBJS) $(CNTLOBJS)
	$(FC) $(FFLAGS) $(OPTFLAGS) $(INCS) $(MICFLAGS) -o $(POP) pop.f90 smooth.o varinfo.o get_tgtgrid.o gaulat.o $(UTILOBJS) $(CNTLOBJS) $(LIBS) $(OMPFLAG_PREPPOST)

$(POPREADINIT): pop_read_init.f90 smooth.o varinfo.o get_tgtgrid.o gaulat.o $(DEPLIBS) $(UTILOBJS) $(CNTLOBJS)
	$(FC) $(FFLAGS) $(OPTFLAGS) $(INCS) $(MICFLAGS) -o $(POPREADINIT) pop_read_init.f90 smooth.o varinfo.o get_tgtgrid.o gaulat.o $(UTILOBJS) $(CNTLOBJS) $(LIBS) $(OMPFLAG_PREPPOST)

$(SLINT_INIT): slint_init_save.f90 get_tgtgrid.o gaulat.o $(DEPLIBS) $(UTILOBJS) $(CNTLOBJS)
	$(FC) $(FFLAGS) $(OPTFLAGS) $(INCS) $(MICFLAGS) $(GPTL_CPPFLAGS) -o $(SLINT_INIT) slint_init_save.f90 get_tgtgrid.o gaulat.o $(UTILOBJS) $(CNTLOBJS) $(LIBS) $(OMPFLAG_PREPPOST) $(GPTL_LDFLAGS)

pop.f90 pop_read_init.f90 slint_init_save.f90: %.f90: %.F90
	$(CPP) $(CPPFLAGS) $< > $@

clean:
	$(RM) *.o *.mod
