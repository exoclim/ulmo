FC=gfortran #fortran compiler
FFLAGS=-O3 -Wall -Wextra -std=f2008 #optimization (for level 3) flags, compiler warnings and the strictest adherence to the latest standards
SRC=NAMELIST.f90 WRITE_READ_DATA.f90 DEGREE_TO_RADIAN.f90 MASS_FLUXES.f90 ULMO.f90
OBJ=${SRC:.f90=.o} #substitute .f90 with .o

%.o: %.f90 #wildcard rule, creation of *.o depends on *.f90
	$(FC) $(FFLAGS) -o $@ -c $<

ULMO: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) 

clean: #cleans all the old compilation files
	@rm -f *.mod *.o ULMO
