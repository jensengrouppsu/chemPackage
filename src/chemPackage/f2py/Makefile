# =============================================================================
# Let's start with the declarations
# =============================================================================

# =============================================================================
# The compiler(s)
LIFX := $(shell which ifx 2>/dev/null)
ifeq ($(LIFX), )
	F90 = gfortran
else
	F90 = ifx
endif

# F90  = gfortran

# Release 
F2PY = FC="${F90}"  python -m numpy.f2py
F90FLAGS =  -heap-arrays -lpthread 
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

ifeq ($(MKLROOT), )
	# No MKL
	F90FLAGS += -lopenblas -lblas -llapack
else
	F90FLAGS += -lirc -limf -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lmkl_rt -lm -ldl
	INC += -I$(MKLROOT)/include
endif

# NOTE: ifort somehow fails sanity chekc of f2py with meson backend 
F2PYFLAGS =  --backend meson --f90flags=${F90FLAGS}

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
	${F2PY} ${F2PYFLAGS} ${INC} -m coordinate_functions \
            -c calc_bonds.f90 calc_dist.f90 \
               minmax_dist.f90 translate.f90 rotate.f90

tests: tests.so

tests.so: find_equal.f90 constants.mod
	${F2PY} ${F2PYFLAGS} ${INC} -m tests -c find_equal.f90

drawingMath: drawingMath.so

drawingMath.so: drawingMath.f90 tensors.f90 constants.mod 
	${F2PY} ${F2PYFLAGS} ${INC} -m drawingMath -c drawingMath.f90 tensors.f90

tensors: tensors.so

tensors.so: tensors.f90 constants.mod
	${F2PY} ${F2PYFLAGS} ${INC} -m tensors -c tensors.f90

constants.mod: Constants.f90
	${F90} ${INC} ${F90FLAGS} -c Constants.f90


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
