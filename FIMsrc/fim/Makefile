# fim Makefile

include ../macros.make

.PHONY: all clean

all:
	$(MAKE) $(GMAKEMINUSJ) -C column FFLAGS="$(COLFLAGS)"
	$(MAKE) $(GMAKEMINUSJ) -C column_chem FFLAGS="$(CHEMFLAGS)"
	$(MAKE) $(GMAKEMINUSJ) -C wrfphys
	(cd horizontal && ./copy.ksh $(HW))
	$(MAKE) $(GMAKEMINUSJ) -C horizontal

clean:
	$(MAKE) -C column      clean
	$(MAKE) -C column_chem clean
	$(MAKE) -C wrfphys     clean
	$(MAKE) -C horizontal  clean
