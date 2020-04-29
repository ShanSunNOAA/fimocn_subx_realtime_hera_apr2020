# cntl Makefile

include ../macros.make

FLAGS   = $(FFLAGS) $(INCS) $(MICFLAGS)
INCS    = -I $(UTILDIR)
LIBCNTL = $(LIBDIR)/libcntl.a
OBJS    = module_control.o units.o module_constants.o module_wrf_control.o fimnamelist.o
UTILDIR = ../utils

.PHONY: all clean

%.o: %.F90
	$(FC) -c $(FLAGS) $<

all: $(LIBCNTL)

$(LIBCNTL): $(OBJS)
	$(AR) ruv $@ $(OBJS)

module_control.o: module_wrf_control.o units.o fimnamelist.o
fimnamelist.o: fimnamelist.exec

clean:
	$(RM) *.o *.mod *.f90
