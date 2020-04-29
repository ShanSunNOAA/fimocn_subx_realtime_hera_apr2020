# gribio Makefile

# This makefile builds library of subroutines that writes grib files.  

include ../../macros.make

FFLAGS	  = $(FFLAGS) $(MICFLAGS)
LIBGRIBIO = $(LIBDIR)/libgribio.a
OBJS      = $(addsuffix .o, $(basename $(SRCS)))
SRCS      = $(wildcard *.F90 *.c)

all: $(LIBGRIBIO)

$(LIBGRIBIO): DEPENDENCIES $(OBJS)
	$(AR) ruv $@ $(OBJS)

%.o: %.F90
	$(FC) -c $(FFLAGS) $(MICFLAGS) $<

%.o: %.c
	$(CC) $(CFLAGS) -c $(GRIBIO_CPPFLAGS) $(MICFLAGS) $<

DEPENDENCIES: $(wildcard *.F90)
	$(RM) Filepath Srcfiles
	echo "." > Filepath
	ls *.F90 > Srcfiles
	$(MKDEPENDS) -m Filepath Srcfiles > $@

-include DEPENDENCIES

clean:
	$(RM) *.o *.mod DEPENDENCIES Filepath Srcfiles
