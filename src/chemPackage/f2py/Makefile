# =============================================================================
# Let's start with the declarations
# =============================================================================

# =============================================================================
# The compiler(s)
F90  = gfortran

# Release 
F2PY = FC=gfortran f2py
F90FLAGS = -FR -fPIC -heap-arrays -fopenmp -pthread -lifcore
# F90OPT   = -ipo -O3 -unroll -inline
# F90ARCH  = -m64
# -ipo      : Interprocedural optimization - for code with more than one file
# -O3       : Optimize for maximum speed, may take longer to compile
# -unroll   : Unroll loops - this removes some overhead associated with looping
# -inline   : Places functions inside main programs
# -FR       : Free format

# debug
#F2PY = f2py --debug
#F90FLAGS = -FR -fPIC -heap-arrays -openmp -threads -lifcore
#F90OPT = -O0 -g -traceback -check all 
#F90ARCH  = -m64


# NOTE: ifort somehow fails sanity chekc of f2py with meson backend 
F2PYFLAGS = --f90flags='${F90FLAGS}' --f77flags='${F90FLAGS}' --backend meson
# --fcompiler=intelem : Use intel's compiler

# =============================================================================
# Libraries
ifeq ($(MKLROOT), )
	MATHLIB = -lblas
else
	MATHLIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core-lpthread
	INC += -I$(MKLROOT)/include
endif

# MKL = -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
INC += -I$(CURDIR)
# =============================================================================

# =============================================================================
# Building rules
# =============================================================================

# "make" builds all
all: coordinate_functions tests drawingMath tensors

coordinate_functions: coordinate_functions.so

coordinate_functions.so: calc_bonds.f90 calc_dist.f90 minmax_dist.f90 \
                         translate.f90 rotate.f90 constants.mod 
	${F2PY} ${F2PYFLAGS} ${INC} ${MATHLIB} -m coordinate_functions \
            -c calc_bonds.f90 calc_dist.f90 \
               minmax_dist.f90 translate.f90 rotate.f90

tests: tests.so

tests.so: find_equal.f90 constants.mod
	${F2PY} ${F2PYFLAGS} ${INC} ${MATHLIB} -m tests -c find_equal.f90

drawingMath: drawingMath.so

drawingMath.so: drawingMath.f90 tensors.f90 constants.mod 
	${F2PY} ${F2PYFLAGS} ${INC} ${MATHLIB} -m drawingMath -c drawingMath.f90 tensors.f90

tensors: tensors.so

tensors.so: tensors.f90 constants.mod
	${F2PY} ${F2PYFLAGS} ${INC} ${MATHLIB} -m tensors -c tensors.f90

constants.mod: Constants.f90
	${F90} ${F90FLAGS} ${INC} ${MATHLIB} -c Constants.f90


# =============================================================================
# Special rules
# =============================================================================

#  Make sure we don't remove this by accident if interrupted at the wrong time.
.PRECIOUS: Makefile

# Utility targets
.PHONY: sources install clean cleanall

# Prints to screen the source files that this program depends on
sources:
	@echo ${SRC}

# Removes objects
clean:
	@rm -vf *.o *.mod

# Removes executable
cleanall: clean
	@rm -vf *.so
