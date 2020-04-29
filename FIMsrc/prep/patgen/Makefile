# Makefile for pattern genertor
include ../../macros.make
#SHELL           =       /bin/sh
.SUFFIXES: .F .i .o .f90 .c

#### Architecture specific settings ####

F90             = ifort
LD              = ifort
#F77             = ifort
#CC              = gcc
#CPP             = cpp -P
#CCFLAGS         = -O3
#CPPFLAGS        =
#NETCDFPATH      = /apps/netcdf/3.6.3-intel

FCFFLAGS        = -assume byterecl -heap-arrays -c
LDFLAGS         = -assume byterecl -heap-arrays
LIB_LOCAL        = libs
LIBPATH         = ../../lib
OPTIMIZATION    = -O3 -xHOST
DEBUG           = -traceback #-C

# Define modules and associated libraries

MODS = kinds.o                                       \
       physcons.o                                    \
       shtns.o                                       \
       patterngenerator.o

PROGOBJ = patgen_icos.o
PROG = $(BINDIR)/patgen_icos

# Define compile time variables

#INCLUDE = -I$(NETCDFPATH)/include -I libs/slint
INCLUDE = -I ../incmod
LIBS   = -L$(LIB_LOCAL)/ -lshtns -L$(LIBPATH) -lw3_4 -lsigio_4 -lslint -mkl=sequential -openmp 
#LIBS   = -L$(LIB_LOCAL)/ -lshtns -lfftw3 -L$(LIBPATH) -lw3_4 -lsigio_4 -lslint -mkl=sequential -openmp 
LIBRARIES = $(LIB_LOCAL)/libw3_4.a $(LIB_LOCAL)/libshtns.a 

# Define executables

all: $(PROG)  

$(PROG): $(LIBRARIES) $(MODS) $(PROGOBJ)  
	$(LD) $(LDFLAGS) $(PROGOBJ) $(INCLUDE) $(MODS) $(LIBS) $(OPTIMIZATION) $(DEBUG) -lm -o $@

$(LIBRARIES): 
	make -C $(LIB_LOCAL)
	
clean:
	rm *.o *.mod patgen_icos
	make -C $(LIB_LOCAL) clean

# Define make rules
.f90.mod:
	$(F90) $(FCFFLAGS) $(INCLUDE) $(OPTIMIZATION) $(DEBUG) -c $*.f90

.f90.o:
	$(F90) $(FCFFLAGS) $(INCLUDE) $(OPTIMIZATION) $(DEBUG) -c $*.f90

.f.o:
	$(F77) $(CPPFLAGS) $(FCFFLAGS) $(INCLUDE) $(OPTIMIZATION) $(DEBUG) -c $*.f

.F.o:
	$(F77) $(CPPFLAGS) $(FCFFLAGS) $(INCLUDE) $(OPTIMIZATION) $(DEBUG) -c $*.F
