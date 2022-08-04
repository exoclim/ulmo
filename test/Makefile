FC=gfortran #fortran compiler
FFLAGS=-O3 -Wall -Wextra -std=f2008 -lgsl -lgslcblas -lm
SRC=./test_fgsl.F90
OBJ=${SRC:/src/.f90=.o} #substitute .f90 with .o

%.o: %.f90 #wildcard rule, creation of *.o depends on *.f90
	$(FC) $(FFLAGS) -o $@ -c $<

test_fgsl: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) 

clean: #cleans all the old compilation files
	@rm -f test.mod test_fgsl.mod *.o test_fgsl
