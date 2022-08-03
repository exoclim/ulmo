FC=gfortran #fortran compiler
FFLAGS=-O3 -Wall -Wextra -std=f2008 #optimization (for level 3) flags, compiler warnings and the strictest adherence to the latest standards
SRC=./src/namelist.F90  ./src/write_read_data.F90 ./src/degree_to_radian.F90 ./src/height_of_slab.F90 ./src/dA_da.F90 ./src/mass_fluxes.F90 ./src/other_fluxes.F90 ./src/div_m.F90 ./src/matrix_calc.F90 ./src/ulmo.F90
OBJ=${SRC:/src/.f90=.o} #substitute .f90 with .o

%.o: %.f90 #wildcard rule, creation of *.o depends on *.f90
	$(FC) $(FFLAGS) -o $@ -c $<

ULMO: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) 

clean: #cleans all the old compilation files
	@rm -f namelist.mod write_read_data.mod degree_to_radian.mod mass_flux.mod da_da.mod div_m.mod height_of_slab.mod other_fluxes.mod matrix_calc.mod *.o ULMO
