# dummygptl Makefile

include ../macros.make

SRCS = $(wildcard *.F90)
OBJS = $(addsuffix .o,$(basename $(SRCS)))

all: libgptl.a

libgptl.a: $(OBJS)
	$(AR) ruv $@ $(OBJS)

%.o: %.F90
	$(FC) -c $(FFLAGS) $(MICFLAGS) $<

clean:
	$(RM) *.o libgptl.a
