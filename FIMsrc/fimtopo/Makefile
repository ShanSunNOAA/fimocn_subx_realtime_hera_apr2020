# fimtopo Makefile

include ../macros.make

CMD      = fimtopo.x
FC       = /opt/lahey/lf6480/bin/f95
FC       = ifort
FFLAGS   =
LD       = /opt/lahey/lf6480/bin/f95
LD       = ifort
LDFLAGS  = -g
LIBS     =	
OBJS     = fimtopo.o

# Lines from here on down should not need to be changed. They are the actual
# rules which make uses to build the executable.

$(CMD): $(OBJS)
	$(LD) -o $(CMD) $(OBJS) $(LIBS) $(LDFLAGS)

clean:
	$(RM) *.o $(CMD) *.mod

%.o: %.f90
	$(info Compiling $<)
	$(FC) -o $@ -c $<
