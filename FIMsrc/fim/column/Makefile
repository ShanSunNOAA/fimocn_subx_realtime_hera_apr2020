# column Makefile

include ../../macros.make

COPYFILES = iw3jdn.f w3fs26.f w3movdat.f w3reddat.f
OBJS      = $(addsuffix .o, $(basename $(SRCS)))
SRCS      = $(wildcard *.f *.F90)

ifeq ($(NEED_SINDCOSD),yes)
  LOCDEF = -DNEED_SINDCOSD
endif

LOCFLAGS = $(FFLAGS)
ifeq ($(DEBUG),yes)
  LOCFLAGS += -g
endif

%.o: %.f
	$(FC) -c -I../../cntl -I../../utils $(GPTL_FFLAGS) $(LOCFLAGS) $(FIXEDFLAG) $(RCWFLAG) $(OPTFLAGS) $(MICFLAGS) $(THREADSTACKFLAG) $<

%.o: %.F90
	$(FC) -c $(GPTL_CPPFLAGS) $(GPTL_FFLAGS) $(LOCDEF) $(LOCFLAGS) $(FREEFLAG) $(RCWFLAG) $(OPTFLAGS) $(MICFLAGS) $(THREADSTACKFLAG) $<

.PHONY: all clean

all: $(OBJS)

$(COPYFILES): %: ../../w3/%
	rsync -a $< $@

-include FIM_COLUMN_DEPENDENCIES

FIM_COLUMN_DEPENDENCIES: $(COPYFILES) $(wildcard *.f *.F90)
	$(RM) Filepath Srcfiles
	echo "." > Filepath
	ls *.f *.F90 > Srcfiles
	$(MKDEPENDS) -m Filepath Srcfiles > $@

clean:
	$(RM) *.o *.mod $(COPYFILES) FIM_COLUMN_DEPENDENCIES
