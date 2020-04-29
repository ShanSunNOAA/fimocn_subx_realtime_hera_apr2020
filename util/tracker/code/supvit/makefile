SHELL=  /bin/sh
ISIZE = 32
RSIZE = 64
#COMP=   xlf
COMP=   ifort
#LIBS=   -L/home/bluevista/dtcrt/CT2007/wrfpostprocV2/lib/ -lw3
LIBS=   -L/tg1/jetmgmt/harrop/WRFV3.0/WPPV3_Intel-10.1/lib -lw3
#LDFLAGS= -bloadmap:loadmap -bnoquiet
LDFLAGS=
#FFLAGS= -O2 -qarch=auto -qintsize=$(ISIZE) -qrealsize=$(RSIZE)
#FFLAGS= -g -O0 -qcheck -qextchk -qarch=auto -qflttrap=ov:zero:inv:enable -qintsize=$(ISIZE) -qrealsize=$(RSIZE)
#FFLAGS= -O2 -qcheck -qarch=auto -qflttrap=ov:zero:inv:enable -qintsize=$(ISIZE) -qrealsize=$(RSIZE)
FFLAGS= -O2 -check all -fpe0 -integer_size $(ISIZE) -real_size $(RSIZE)



supvit:     supvit_main.f supvit_modules.o
	@echo " "
	@echo "  Compiling program that sorts and updates vitals records...."
	$(COMP) $(FFLAGS) $(LDFLAGS) supvit_modules.o supvit_main.f $(LIBS) -o supvit
	@echo " "

supvit_modules.o:   supvit_modules.f
	@echo " "
	@echo "  Compiling the modules....."
	$(COMP) -c supvit_modules.f -o supvit_modules.o
	@echo " "

