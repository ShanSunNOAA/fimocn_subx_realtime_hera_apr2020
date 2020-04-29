# horizontal Makefile

include ../../macros.make
include FIM_HORIZONTAL_OBJS

AR               = ar  # TODO: move into top-level makefile or macros.make
ARFLAGS          = ruv # TODO: move into top-level makefile or macros.make
COLMOBJS         = $(wildcard ../column/*.o)
COLUMN_CHEM_OBJS = $(wildcard ../column_chem/*.o)
DEPLIBS          = $(foreach x,bacio_4 sfcio_4 sigio_4 slint sp_4 sysshare w3_4 ip, $(LIBDIR)/lib$(x).a) 
DOTSMS           = machine.sms namelist_def.sms
FIMEXE           = $(BINDIR)/fim
FIMLIB           = $(LIBDIR)/lib$(FIMLIBBASE).a
FIMLIBBASE       = fim
ICOSIO           = $(wildcard ../../icosio/*.o)
INCS1           = $(addprefix -I../,column column_chem wrfphys)
INCS2           = $(addprefix -I../../,icosio post/pop prep/sfcio prep/sigio prep/incmod w3)
INCUTILSDIR      = -I$(UTILSDIR)
LIBS             = -l$(FIMLIBBASE) -L$(LIBDIR) $(addprefix -l,bacio_4 gribio sfcio_4 sigio_4 slint sp_4 sysshare w3_4 ip) 
LOCFLAGS         = $(FFLAGS) $(FREEFLAG) $(INCS1) $(INCS2) $(INCUTILSDIR)
POSTOBJS         = $(foreach x,varinfo post smooth get_tgtgrid gaulat,../../post/pop/$(x).o)
PPP              = $(SMS)/bin/ppp -I $(INCMPI)
SMSEXTENSION     = _sms
UTILS            = $(UTILSDIR)/headers.o
UTILSDIR         = ../../utils
WRFPHYS_OBJS     = $(filter-out %module_cu_g3.o,$(wildcard ../wrfphys/*.o))
#PPPINCS          = -I ../../prep/lsfit/par
#MKLIB = -lmkl_intel_lp64 -lmkl_core -lmkl_sequential

# Add SMS include directory to search path for parallel builds.

ifeq ($(P),p)
  LOCFLAGS += -I$(SMS)/include
endif

#ifeq ($(P),p)
#  LOCFLAGS += -I ../../prep/lsfit/par
#  DEPLIBS += ../../prep/lsfit/par/liblsfit.a
#  LIBS +=  -L../../prep/lsfit/par/ -llsfit $(MKLIB)
#endif

#ifeq ($(P),s)
#  LOCFLAGS += -I ../../prep/lsfit/ser
#  DEPLIBS += ../../prep/lsfit/ser/liblsfit.a
#  LIBS +=  -L../../prep/lsfit/ser/ -llsfit $(MKLIB)
#endif

# Add -g for debugging.

ifeq ($(DEBUG),yes)
 LOCFLAGS += -g
endif

# Add SERIAL token to CPPFLAGS for serial builds.

ifeq ($(P),s)
  CPPFLAGS += -DSERIAL
endif

# If degree-based trig intrinsics don't exist, set ifdef to define them as stmt
# functions.

ifeq ($(NEED_SINDCOSD),yes)
  CPPFLAGS += -DNEED_SINDCOSD
endif

# Are we going to attach gdb or ddd to a running process?

ifeq ($(ATTACH_DEBUGGER),yes)
  LOCFLAGS += -DATTACH_DEBUGGER
endif

ifeq ($(HW),gpu)
  ifeq ($(P),p)
    LIBSMS = -lsmsgpu
  endif
  LIB_CUDA = $(LIBCUDA)
  OBJS_BASE = $(OBJS_CPU)
  OBJS_XTRA = $(OBJS_GPU)
else
  ifeq ($(OPENMP),yes)
    ifeq ($(P),p)
      LIBSMS = -lsmsomp
    endif
  else
    ifeq ($(P),p)
      LIBSMS = -lsms
    endif
  endif
  LIB_CUDA =
  OBJS_BASE = $(OBJS_CPU) $(OBJS_GPU)
  OBJS_XTRA = 
endif

%.o: %.F90
	$(CPP) $(GPTL_CPPFLAGS) $(CPPFLAGS) $< > $*.f90
ifeq ($(P),s)
	$(FC) -c $(OMPFLAG) $(GPTL_FFLAGS) $(LOCFLAGS) $(OPTFLAGS) $(INCNETCDF) $*.f90 -o $*.o
else
	$(PPP) $(PPPINCS) $*.f90 $*$(SMSEXTENSION).f90
	$(FC) -c $(OMPFLAG) $(GPTL_FFLAGS) $(LOCFLAGS) $(OPTFLAGS) $(INCNETCDF) $*$(SMSEXTENSION).f90 -o $*.o
endif

%.o: %.c
	$(CC) -c -g $(CFLAGS) $(GCCOMPFLAG) $(FORTRAN_UNDERSCORING) $<

.PHONY: all clean gpu

all: gpu $(FIMEXE)

gpu:
ifeq ($(HW),gpu)
	make -f Makefile.serial.gpu 
endif

physics.o: $(DEPLIBS)

fimnamelist.f90: $(addprefix fimnamelist,.F90 .exec)

#JR May want to compile ss2icos.F90 with threading enabled. So need separate rule
ss2icos.o:
	$(CPP) $(GPTL_CPPFLAGS) $(CPPFLAGS) $< > $*.f90
ifeq ($(P),s)
	$(FC) -c $(OMPFLAG) $(OMPFLAG_SS2ICOS) $(GPTL_FFLAGS) $(LOCFLAGS) $(OPTFLAGS) $*.f90 -o $*.o
else
	$(PPP) $*.f90 $*$(SMSEXTENSION).f90
	$(FC) -c $(OMPFLAG) $(OMPFLAG_SS2ICOS) $(GPTL_FFLAGS) $(LOCFLAGS) $(OPTFLAGS) $*$(SMSEXTENSION).f90 -o $*.o
endif

# Conditionally generate .sms module files for sources not processed by ppp.

ifeq ($(P),p)
gfsphys_nl.o: $(DOTSMS)
namelist_def.sms: machine.sms
$(DOTSMS): %.sms: ../column/%.f
	$(PPP) fixed modinfo ../column/$*.f .
endif

# Include (and, if necessary, generate) dependencies makefile.

-include FIM_HORIZONTAL_DEPENDENCIES

FIM_HORIZONTAL_DEPENDENCIES: $(wildcard *.F90)
	$(RM) Filepath Srcfiles
	echo "." > Filepath
	ls *.F90 > Srcfiles
	$(MKDEPENDS) -m -d sms__decomp.o Filepath Srcfiles > $@

# Build the FIM executable.

include FIM_HORIZONTAL_OBJS_TOP

$(FIMEXE): $(FIMLIB) $(OBJS_TOP) $(DEPLIBS)
	$(FC) -o $(FIMEXE) $(OMPFLAG) $(OMPFLAG_SS2ICOS) $(FFLAGS) -L$(SMS)/lib $(LINKFLAGS) $(OBJS_TOP) $(OBJS_XTRA) \
          $(LIBS) $(LIBSMS) $(LINKLIBS) $(GPTL_LDFLAGS) $(LDFLAGS) $(LIBNETCDF) $(LIB_CUDA)

# Build a library of all FIM object files except $(OBJS_TOP).  This is needed
# when building FIM from within the NEMS repository but would not be needed by
# other builds.  However, always building $(FIMLIB) significantly simplifies
# build automation (for example, no additional $(FTNMPI) values need to be
# supported by makefim, special logic is not needed in fim/horizontal/copy.ksh
# to avoid building $(OBJS_TOP) files, additional logic is not needed in the
# top-level Makefile, etc.).  And, since we will always build it, we should
# use it to build $(FIMEXE) to avoid using a separate but equivalent approach
# when building from within NEMS.

$(FIMLIB): $(OBJS_BASE) $(COLMOBJS) $(POSTOBJS) $(COLUMN_CHEM_OBJS) $(WRFPHYS_OBJS) $(UTILS) $(ICOSIO)
	$(AR) $(ARFLAGS) $@ $?
