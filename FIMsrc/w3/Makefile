# w3 Makefile

###############################################################
#
#   AUTHOR:    Gilbert - W/NP11
#
#   DATE:      01/11/1999
#
#   DATE:      03/25/2010
#   Modified by Rosinski (NOAA/ESRL) to use only a Makefile,
#   eliminate the need for a script
#
###############################################################

#
#     Generate a list of object files that corresponds to the
#     list of Fortran ( .f, .F90 ) files in the current directory
#

include ../macros.make

LIB  = $(LIBDIR)/libw3_4.a
OBJS = $(addsuffix .o, $(basename $(SRCS)))
SRCS = $(wildcard *.f *.F *.F90)

.PRECIOUS: %.f90

.PHONY: all clean

%.o: %.F90
	$(FC) $(FFLAGS_NO_DEBUG) $(FREEFLAG) $(MICFLAGS) -c $<

%.o: %F
	$(FC) $(FFLAGS_NO_DEBUG) $(FREEFLAG) $(MICFLAGS) -c $<

%.o: %.f
	$(FC) -c $(FFLAGS) $(FIXEDFLAG) $(MICFLAGS) $<

all: $(LIB) $(BINDIR)/jdate

$(LIB): $(OBJS)
	$(AR) ruv $@ $(OBJS)

$(BINDIR)/jdate: jdate.f90 $(LIB)
	$(FC) $(FFLAGS) $(FREEFLAG) $(MICFLAGS) -o $@ $< -L$(LIBDIR) $(LIB)

errexit.f: errexit.F
	$(CPP) $(CPPFLAGS) $< > $@

errexit.o: errexit.f
	$(FC) -c $(FFLAGS) $(MICFLAGS) $<

jdate.f90: jdate.F90
	$(CPP) $(CPPFLAGS) $< > $@

jdate.o: jdate.f90
	$(FC) -c $(FFLAGS) $(FREEFLAG) $(MICFLAGS) $<

clean:
	$(RM) *.o *.mod
