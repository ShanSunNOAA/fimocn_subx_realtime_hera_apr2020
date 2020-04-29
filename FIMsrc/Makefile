# Makefile

include macros.make

ifeq ($(P),s)
  CPPFLAGS += -DSERIAL
endif

.PHONY: all build clean cleanbelow info

all: info build

info:
	$(info Inside top level Makefile)
	$(info BACIO_FFLAGS    = $(BACIO_FFLAGS))
	$(info BINDIR          = $(BINDIR))
	$(info CC              = $(CC))
	$(info CFLAGS          = $(CFLAGS))
	$(info CHEM_FFLAGS     = $(CHEM_FFLAGS))
	$(info COL_FFLAGS      = $(COL_FFLAGS))
	$(info CPPFLAGS        = $(CPPFLAGS))
	$(info CURDIR          = $(CURDIR))
	$(info FFLAGS_NO_DEBUG = $(FFLAGS_NO_DEBUG))
	$(info FIM_FFLAGS      = $(FIM_FFLAGS))
	$(info FIXEDFLAG       = $(FIXEDFLAG))
	$(info FREEFLAG        = $(FREEFLAG))
	$(info GPTL_FLAGS      = $(GPTL_FLAGS))
	$(info GRIBIO_CPPFLAGS = $(GRIBIO_CPPFLAGS))
	$(info ICOSIO_FFLAGS   = $(ICOSIO_FFLAGS))
	$(info INCNETCDF       = $(INCNETCDF))
	$(info LIBDIR          = $(LIBDIR))
	$(info LIBNETCDF       = $(LIBNETCDF))
	$(info LINKFLAGS       = $(LINKFLAGS))
	$(info LINKLIBS        = $(LINKLIBS))
	$(info MKDEPENDS       = $(MKDEPENDS))
	$(info NEMS            = $(NEMS))
	$(info POST_FFLAGS     = $(POST_FFLAGS))
	$(info PREP_FFLAGS     = $(PREP_FFLAGS))
	$(info SP_FFLAGS       = $(SP_FFLAGS))
	$(info W3_FFLAGS       = $(W3_FFLAGS))

build:
	$(MAKE) FFLAGS="$(GPTL_FLAGS)"                                                     -C dummygptl
	$(MAKE) CPPFLAGS="$(CPPFLAGS)" FFLAGS="$(ICOSIO_FFLAGS) $(OPTFLAGS) $(MICFLAGS)"   -C icosio
	$(MAKE) FFLAGS="$(PREP_FFLAGS)"                                                    -C sys_share
	$(MAKE) $(GMAKEMINUSJ) FFLAGS="$(PREP_FFLAGS)"                                     -C utils
	$(MAKE) FFLAGS="$(PREP_FFLAGS)"                                                    -C cntl
	$(MAKE) $(GMAKEMINUSJ) COPT="$(CFLAGS) $(MICFLAGS)" FOPT="$(BACIO_FFLAGS)"         -C bacio
	$(MAKE) $(GMAKEMINUSJ) FFLAGS="$(W3_FFLAGS)"                                       -C w3
	$(MAKE) FFLAGS="$(PREP_FFLAGS)"                                                    -C prep
	$(MAKE) FFLAGS="$(POST_FFLAGS)"                                                    -C post
	$(MAKE) CHEMFLAGS="$(CHEM_FFLAGS)" COLFLAGS="$(COL_FFLAGS)" FFLAGS="$(FIM_FFLAGS)" -C fim

clean:
	$(RM) -r ../FIMsrc_*

# Do not use this unless you know what you are doing.

cleanbelow:
	touch macros.make
	$(RM) bin/*
	$(RM) lib/*
	$(MAKE) -C bacio  clean
	$(MAKE) -C cntl   clean
	$(MAKE) -C fim    clean
	$(MAKE) -C icosio clean
	$(MAKE) -C post   clean
	$(MAKE) -C prep   clean
	$(MAKE) -C utils  clean
	$(MAKE) -C w3     clean
	$(RM) macros.make
