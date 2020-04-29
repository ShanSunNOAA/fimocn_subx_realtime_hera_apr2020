# slint Makefile

SHELL = /bin/sh

include ../../macros.make

FLAGS       = $(FFLAGS)
LIBSYSSHARE = $(LIBDIR)/libsysshare.a
SLINTEXE    = $(BINDIR)/slint
SLINTLIB    = $(LIBDIR)/libslint.a
INCMOD	 = ../incmod
#JR Add OPENMP flag for slint to be thread-safe
#JR TODO: construct things so the flag doesn't need to appear TWICE!
FLAGS = $(FFLAGS) -I$(INCMOD) $(OMPFLAG_SS2ICOS) $(OMPFLAG) $(MICFLAGS) $(OPTFLAGS)

#JR: Blank out all suffix rules. Mac does some weird .s stuff without it.

.SUFFIXES:
.SUFFIXES: .F90 .o

.F90.o:
	$(FC) $(FLAGS) -c $<

KDOBJ = kd.o
SLINTOBJ = slint.o uvxyz.o uv2uv.o

all: $(SLINTLIB) $(SLINTEXE)

$(SLINTLIB): DEPENDENCIES $(KDOBJ) $(SLINTOBJ) 
	$(AR) ruv $(SLINTLIB) $(SLINTOBJ) $(KDOBJ)
	mv -f *.mod $(INCMOD)

$(SLINTEXE): DEPENDENCIES slintest.o $(SLINTLIB)
	$(FC) $(FLAGS) $(MICFLAGS) -o $(SLINTEXE) slintest.o $(SLINTLIB) $(LIBSYSSHARE)

DEPENDENCIES:
	$(RM) -f Filepath Srcfiles
	echo "." > Filepath
	ls -1 *.F90 > Srcfiles
	$(MKDEPENDS) -m Filepath Srcfiles > $@

-include DEPENDENCIES

clean:
	$(RM) *.o *.mod *_cpp.F90 DEPENDENCIES Filepath Srcfiles

# Ning Wang, Feb 2006
