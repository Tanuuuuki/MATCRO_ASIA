PROGRAM = MATCRO_Asia
FC      = ifort
FILES   = MATCRO.f90 \
	SUB_IO.f90 \
	SUB_CROP.f90 \
	SUB_PHSYN.f90 \
	SUB_SOIL.f90 \
	SUB_RAD.f90 \
	SUB_TINTERP.f90 \
	FNC.f90
OBJS    = $(FILES:%.f90=%.o)
FNCDF   = /home/appl/netcdf-fortran-4.5.4
CNCDF 	= /home/appl/netcdf-c-4.8.1

%.o: %.f90
	$(FC) -o $@ -CB -c -fpe0 -traceback -g -I$(FNCDF)/include -I$(CNCDF)/include $<

$(PROGRAM): $(OBJS)
	$(FC) -o $@ -CB -fpe0 $^ -L$(FNCDF)/lib -L$(CNCDF)/lib -lnetcdf -lnetcdff

clean:
	rm -f *.o

.SUFFIXES:
.SUFFIXES:.o .f90
