# utils Makefile

include ../macros.make

BINS        = extract_atcf get_num_cores reduce
CNTLDIR    = ../cntl
CNTLSRCS   = $(addprefix $(CNTLDIR)/,units.F90 fimnamelist.F90 fimnamelist.exec)
COMMON      = module_initial_chem_namelists.o fimnamelist.o units.o
INSTALLS    = $(addprefix $(BINDIR)/,$(BINS))
LIBSYSSHARE = $(LIBDIR)/libsysshare.a
MKBIN       = $(FC) $(FFLAGS) $(MICFLAGS) -o $@ $(addsuffix .o,$(notdir $@))
SYNCS       = $(CNTLSRCS)

.PRECIOUS: %.f90

%.f90: %.F90
	$(CPP) $(CPPFLAGS) $< > $@

%.o: %.f90
	$(FC) -c $(DEBUG_FLAGS) $(FFLAGS) $(FREEFLAG) $(MICFLAGS) $<

all: $(INSTALLS) headers.o module_libmf.o

$(INSTALLS): $(BINS)
	install $(notdir $@) $(BINDIR)

extract_atcf reduce: %: %.o
	$(MKBIN)

get_num_cores: %: %.o
	$(MKBIN) wtinfo.o $(COMMON) $(LIBSYSSHARE)

wtinfo.o: wtinfo.f90
	$(FC) -c $(DEBUG_FLAGS) $(FFLAGS) $(FREEFLAG) $(LIBSYSSHARE) $<

fimnamelist.f90: units.o $(addprefix fimnamelist,.F90 .exec)

DEPENDENCIES: $(SYNCS)
	rsync -av $(SYNCS) .
	$(RM) Filepath Srcfiles
	echo "." > Filepath
	ls *.F90 > Srcfiles
	$(MKDEPENDS) -m Filepath Srcfiles > $@

-include DEPENDENCIES

clean:
	$(RM) *.o *.mod DEPENDENCIES Filepath Srcfiles $(SYNCS)
