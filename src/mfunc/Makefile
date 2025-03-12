#
# =============================================================================
# Let's start with the declarations
# =============================================================================

# =============================================================================
# The compiler(s)
F90 = gfortran
# =============================================================================
# =============================================================================
# Flags
F2PY = FC=gfortran f2py
F90FLAGS = -FR -fPIC -heap-arrays -fopenmp -pthread -lifcore
# F90FLAGS += -FR -fPIC
# F90OPT += -ipo -O3 -unroll -inline
# F90ARCH += -m64
# -ipo      : Interprocedural optimization - for code with more than one file
# -O3       : Optimize for maximum speed, may take longer to compile
# -unroll   : Unroll loops - this removes some overhead associated with looping
# -inline   : Places functions inside main programs
# -FR       : Free format
# NOTE: ifort somehow fails sanity chekc of f2py with meson backend 
F2PYFLAGS = --f90flags='${F90FLAGS}' --f77flags='${F90FLAGS}' --backend meson
# --fcompiler=intelem : Use intel's compiler
#  f2py --fcompiler=intelem --f90flags=' -FR -fPIC'  --arch='-m64' --opt='-ipo -O3 -unroll -inline'
# =============================================================================
# Libraries
# MKLROOT = /usr/global/intel/Compiler/2011.11/339/mkl
#MKL = -L$(MKLROOT)/lib/intel64/ -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm
ifeq ($(MKLROOT), )
	MATHLIB = -lblas
else
	MATHLIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
	INC += -I$(MKLROOT)/include
endif

# MKL = -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
INC += -I$(CURDIR)
# =============================================================================

# =============================================================================
# Building rules
# =============================================================================

# "make" builds all
all: erf mvmult

# General rule for building executable from objects. 
# $@ is the name of the target (in this case the executable)
mvmult: mvmult.so

mvmult.so: mvmult.f90 constants.mod
	${F2PY} ${F2PYFLAGS} ${INC} ${MATHLIB} -m mvmult -c mvmult.f90

erf: erf.so

erf.so: erf.f90 constants.mod
	${F2PY} ${F2PYFLAGS} ${INC} ${MATHLIB} -m erf -c erf.f90

constants.mod: Constants.f90
	${F90} ${F90FLAGS} ${INC} ${MATHLIB} -c Constants.f90

# =============================================================================
# Special rules
# =============================================================================

#  Make sure we don't remove this by accident if interrupted at the wrong time.
.PRECIOUS: Makefile

# Utility targets
.PHONY: sources clean cleanall

# Prints to screen the source files that this program depends on
sources:
	@echo ${SRC}

# Removes objects
clean:
	@rm -vf *.o *.mod

# Removes executable
cleanall: clean
	@rm -vf *.so
