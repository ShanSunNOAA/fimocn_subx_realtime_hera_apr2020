# prep Makefile

include ../macros.make

all:
	$(MAKE) -C sfcio
	$(MAKE) -C sigio $(GMAKEMINUSJ) FFLAGS="$(SP_FFLAGS)"
	$(MAKE) -C slint $(GMAKEMINUSJ)
	$(MAKE) -C lsfit $(GMAKEMINUSJ)
	$(MAKE) -C grid  $(GMAKEMINUSJ)
# ifneq ($(DOPATGEN),no)
# 	$(MAKE) -C patgen 
# endif
	$(MAKE) -C sp    $(GMAKEMINUSJ) FFLAGS="$(SP_FFLAGS)"

ifneq ($(DOSSFC2ICOS),no)
	$(MAKE) -C ssfc2icos
endif
	$(MAKE) -C gfsenkf

clean:
	$(MAKE) -C gfsenkf clean
	$(MAKE) -C grid clean
	$(MAKE) -C sfcio clean
	$(MAKE) -C sigio clean
	$(MAKE) -C slint clean
	$(MAKE) -C lsfit clean
# 	$(MAKE) -C patgen clean
	$(MAKE) -C sp clean
	$(MAKE) -C ssfc2icos clean
