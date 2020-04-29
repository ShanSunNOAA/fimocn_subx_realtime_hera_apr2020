lib=baopenr.ansi.standard.o  \
	baread.ansi.standard.o  \
	baopenw.ansi.standard.o \
        wryte.o

balib.a:	$(lib)
	touch balib.a
	rm balib.a
	ar -rc balib.a $(lib)
	rm *.o

.f.o:
	$(FC) -c $(FFLAGS) $<
