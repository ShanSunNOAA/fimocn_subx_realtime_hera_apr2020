# ssfc2icos Makefile

include ../../macros.make

DEPLIBS    = $(LIBW3) $(LIBSP) $(LIBSIGIO) $(LIBSFCIO) $(LIBSLINT) $(LIBSYSSHARE)
FLAG1      = $(FFLAGS)
FLAGS      = $(FLAG1) $(INCS) $(DEBUG_FLAGS) $(MICFLAGS)
INCS       = -I ../incmod -I../../utils -I../../cntl -I../sfcio $(INCNETCDF)
LIBS       = -L$(LIBDIR) -lw3_4 -lsp_4 -lsigio_4 -lsfcio_4 -lslint
LIBSFCIO   = $(LIBDIR)/libsfcio_4.a
LIBSIGIO   = $(LIBDIR)/libsigio_4.a
LIBSLINT   = $(LIBDIR)/libslint.a
LIBSP      = $(LIBDIR)/libsp_4.a
LIBSYSSHARE= $(LIBDIR)/libsysshare.a
LIBW3      = $(LIBDIR)/libw3_4.a
NEWNAME    = $(BINDIR)/newname.exe
SSFC2ICOS  = $(BINDIR)/ssfc2icos

#mf 20130206 -- modules needed by mkmodisland

UTILDIR    = ../../utils
CNTLDIR    = ../../cntl

# not really needed outside fime
#GBLBDIR    = ../../fim/horizontal

UTILOBJS   = $(UTILDIR)/headers.o \
	     $(UTILDIR)/module_initial_chem_namelists.o \
	     $(UTILDIR)/module_libmf.o

CNTLOBJS   = $(CNTLDIR)/module_control.o \
	     $(CNTLDIR)/units.o \
	     $(UTILDIR)/fimnamelist.o \
	     $(CNTLDIR)/module_wrf_control.o \
	     $(CNTLDIR)/module_constants.o

#GBLBOBJS   = $(GBLBDIR)/global_bounds.o

NEWOBJS    = $(UTILOBJS) $(CNTLOBJS)

SSFCLIBS   = $(LIBS) $(LIBSLINT) $(LIBSYSSHARE) $(LIBNETCDF)

%.o: %.F90
	$(FC) -c $(FLAGS) $<

all: $(SSFC2ICOS) $(NEWNAME)

$(NEWNAME): newname.F90 $(DEPLIBS) $(UTILDIR)/headers.o $(CNTLDIR)/fimnamelist.o
	$(FC) $(FLAGS) -o $(NEWNAME) $(INCS) -I../../cntl/incmod $(GPTL_CPPFLAGS) newname.F90 $(LIBS) -lcntl $(UTILOBJS) $(OMPFLAG) $(OMPFLAG_SS2ICOS) $(GPTL_LDFLAGS) $(LIBSYSSHARE)

#MF -- don't need global_bounds.o because not in fim; new .o must be before $SFCLIBS
$(SSFC2ICOS): ssfc2icos.F90 $(DEPLIBS) read_mtnvar.o mkmodisland.o $(NEWOBJS)
	$(FC) $(FLAGS) -o $(SSFC2ICOS) $(INCS) $(GPTL_CPPFLAGS) ssfc2icos.F90  read_mtnvar.o mkmodisland.o $(SSFCLIBS) $(NEWOBJS) $(LIBSYSSHARE) $(OMPFLAG_PREPPOST) $(GPTL_LDFLAGS)

read_mtnvar.o: read_mtnvar.F90
	$(FC) -c $(FLAGS) read_mtnvar.F90

mkmodisland.o: mkmodisland.F90
	$(FC) -c $(FLAGS) $(INCNETCDF) $(LIBNETCDF) mkmodisland.F90

clean:
	$(RM) *.o *.mod
