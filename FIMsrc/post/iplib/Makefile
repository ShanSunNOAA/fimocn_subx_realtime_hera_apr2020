# iplib Makefile
 
include ../../macros.make

IPLIBOBJ =   \
  gdswiz.o   \
  gdswiz00.o \
  gdswiz01.o \
  gdswiz03.o \
  gdswiz04.o \
  gdswiz05.o \
  gdswizc9.o \
  gdswizca.o \
  gdswizcb.o \
  gdswzd.o   \
  gdswzd00.o \
  gdswzd01.o \
  gdswzd03.o \
  gdswzd04.o \
  gdswzd05.o \
  gdswzdc9.o \
  gdswzdca.o \
  gdswzdcb.o

%.o: %.F90
	$(FC) -c -g $(FFLAGS) $(OPTFLAGS) $(INCS) $<

%.o: %.f
	$(FC) -c $(FFLAGS) $<

$(LIBDIR)/libip.a: $(IPLIBOBJ)
	$(AR) ruv $@ $^
