# sfcio Makefile

include ../../macros.make

#JR sfcio_module takes a long time to compile and does tiny computation so use -O0
ifdef NAG
FLAGS    = $(FFLAGS) $(DEBUG_FLAGS) $(GPTL_FFLAGS) $(MICFLAGS) -O0
else
FLAGS    = $(FFLAGS_NO_DEBUG) $(DEBUG_FLAGS) $(GPTL_FFLAGS) $(MICFLAGS) -O0
endif
LIBSFCIO = $(LIBDIR)/libsfcio_4.a

.PHONY: all clean

all: $(LIBSFCIO)

$(LIBSFCIO): $(LIBSFCIO)(sfcio_module.o)

%.o: %.F90
	$(FC) -c $(FLAGS) $<

clean:
	$(RM) *.o *.mod
