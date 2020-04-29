################################################################################
# This makefile was produced by fmjif at 17:21:11 on 97/02/24.
SHELL=	/bin/sh
CMD=	copygb
# FOPTS=	-O2 -i4 -r8
FOPTS=	-O -m64 -i4 -r8
LOPTS=
INCS=	
OBJS=	copygb.o
#LIBS= -L/nwprod/lib -lw3_d -lbacio_4 -lip_d -lsp_d	
# LIBS= -L../libs -lbacio_8 -lip_8 -lsp_8 -L/home/tpm/tracker/libs -lw3_d
LIBS= -L/mnt/tg1/projects/wrf-dte/tracker/code/copygb_fd/copygb_gfdl/  -lip -lsp -lw3
################################################################################
$(CMD):	$(OBJS)
	ifort $(LOPTS) $(OBJS) $(LIBS) -o $(CMD)
$(OBJS):	$(INCS)
copygb.o:	copygb.f
	ifort -c $(FOPTS) $*.f
clean:	
	-rm -f $(OBJS)
