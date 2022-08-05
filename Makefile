.PHONY: all clean

PREFIX=${CONDA_PREFIX}
FC=gfortran #fortran compiler
FFLAGS=-O3 -Wall -Wextra -std=f2008 
INCLUDE=-I${PREFIX}/include/fgsl -I${PREFIX}/include
LIBS=-L${PREFIX}/lib -lfgsl -lgsl -lgslcblas -lm -lcblas -lm

EXE=ulmo.out
SRC=./src/namelist.F90 ./src/read_data.F90 ./src/degree_to_radian.F90 ./src/height_of_slab.F90 ./src/dA_da.F90 ./src/mass_fluxes.F90 ./src/other_fluxes.F90 ./src/div_m.F90 ./src/mat_test.F90 ./src/matrix_calc.F90 ./src/process_output_data.F90 ./src/ulmo.F90
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
