SHELL = /bin/sh

include ../../macros.make

FLAGS       = $(FFLAGS)
LIBLSFIT    = ./liblsfit.a
FLAGS = $(FFLAGS) $(MICFLAGS)

.SUFFIXES:
.SUFFIXES: .F90 .o

%.o: %.F90
	$(FC) -c $(FLAGS) -o $@ $<
%.o: %.f
	$(FC) -c $(FLAGS) -o $@ $<

SVDOBJ = la_aux.o svd.o
BLASOBJ = blas.o
LSFITOBJ = LSfitting.o

all: $(LIBLSFIT)

$(LIBLSFIT): $(SVDOBJ) $(BLASOBJ) $(LSFITOBJ)
	$(AR) ruv $(LIBLSFIT) $(SVDOBJ) $(BLASOBJ) $(LSFITOBJ) 

clean:
	$(RM) *.o *.mod *.a

