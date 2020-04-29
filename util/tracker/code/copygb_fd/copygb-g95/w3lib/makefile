w3lib=	errexit.o\
	errmsg.o\
	getgbem.o\
	getgbemh.o\
	getgb1r.o\
	putgben.o\
	w3tagb.o\
	getgi.o\
	getgir.o\
	getbit.o\
	getgb1s.o\
	getgbem.o\
	skgb.o\
	ixgb.o\
	fparser.o fparsei.o\
	instrument.o lengds.o \
	w3utcdat.o w3difdat.o w3reddat.o\
	iw3jdn.o\
	pdsens.o\
	pdseup.o\
        gbytes.o gbytes_char.o\
	mova2i.o\
	xmovex.o\
	r63w72.o\
	w3locdat.o\
	w3doxdat.o\
	w3fs26.o w3fi01.o w3fi58.o w3fi59.o w3fi63.o w3fi68.o w3fi71.o w3fi72.o w3fi73.o w3fi74.o w3fi75.o w3fi76.o w3fi82.o w3fi83.o




w3lib.a:	$(w3lib)
	touch w3lib.a
	rm w3lib.a
	ar -rc w3lib.a $(w3lib)

.f.o:
	$(FC) -c $(FFLAGS) $<

.c.o:
	$(CC) -c $(CFLAGS) $<
