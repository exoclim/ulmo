.PHONY: all clean

PREFIX=${CONDA_PREFIX}
FC=gfortran #fortran compiler
FFLAGS=-O3 -Wall -Wextra -std=f2008 -fcheck=all -g -fall-intrinsics  
INCLUDE=-I${PREFIX}/include/fgsl -I${PREFIX}/include
LIBS=-L${PREFIX}/lib -lfgsl -lgsl -lgslcblas -lm -lcblas -lm

<<<<<<< HEAD
EXE=ulmo.out
#Uncomment 1st SRC if want fgsl solver
#SRC=./src/namelist.F90 ./src/read_data.F90 ./src/degree_to_radian.F90 ./src/height_of_slab.F90 ./src/mass_fluxes.F90 ./src/heat_fluxes.F90 ./src/dA_da.F90 ./src/div_m.F90 ./src/matrix_calc.F90 ./sr#c/memory_usage.F90 ./src/process_output_data.F90 ./src/ulmo_fgsl.F90
SRC=./src/namelist.F90 ./src/read_data.F90 ./src/degree_to_radian.F90 ./src/height_of_slab.F90 ./src/mass_fluxes.F90 ./src/heat_fluxes.F90 ./src/dA_da.F90 ./src/div_m.F90 ./src/process_output_data.F90 ./src/new_temp_calc_fe.F90 ./src/ulmo_fe.F90

=======
EXE=main.out
SRC=./src/constants.F90 ./src/read_data.F90 ./src/degree_to_radian.F90 ./src/height_of_slab.F90 ./src/heat_fluxes.F90 ./src/matrix_calc.F90 ./src/process_output_data.F90 ./src/main.F90
>>>>>>> 54079431224abd92ef840aca5829df8ea309df91
OBJ=${SRC:/src/.f90=.o} #substitute .f90 with .o

help:
	@echo 'Makefile for the test code                                     '
	@echo '                                                               '
	@echo 'Usage:                                                         '
	@echo '    make all                                   Compile the code'
	@echo '    make run                                 Run the executable'
	@echo '    make clean           Remove the executable and object files'
	@echo '                                                               '

all: $(EXE)

run: ${EXE}
	./${EXE}

clean: # cleans all the old compilation files
	@rm -f *.mod *.o *.out

%.o: %.f90 # wildcard rule, creation of *.o depends on *.f90
	$(FC) $(FFLAGS) -o $@ -c $< ${INCLUDE}

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) ${INCLUDE} ${LIBS}
