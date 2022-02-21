#
# =============================================================================
# Let's start with the declarations
# =============================================================================

# =============================================================================
# The compiler(s)
F90 = ifort
F2PY = f2py
# =============================================================================
# =============================================================================
# Flags
F90FLAGS += -FR -fPIC
F90OPT += -ipo -O3 -unroll -inline
F90ARCH += -m64
# -ipo      : Interprocedural optimization - for code with more than one file
# -O3       : Optimize for maximum speed, may take longer to compile
# -unroll   : Unroll loops - this removes some overhead associated with looping
# -inline   : Places functions inside main programs
# -FR       : Free format
F2PYFLAGS += --fcompiler=intelem --f90flags='${F90FLAGS}' --arch='${F90ARCH}' \
             --opt='${F90OPT}'
# --fcompiler=intelem : Use intel's compiler
#  f2py --fcompiler=intelem --f90flags=' -FR -fPIC'  --arch='-m64' --opt='-ipo -O3 -unroll -inline'
# =============================================================================
# Libraries
MKLROOT = /usr/global/intel/Compiler/2011.11/339/mkl
#MKL = -L$(MKLROOT)/lib/intel64/ -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm
MKL = -L$(MKLROOT)/lib/intel64/ -lmkl_rt -lpthread -lm
INC += -I$(CURDIR) -I$(MKLROOT)/include
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
	${F2PY} ${F2PYFLAGS} ${INC} ${MKL} -m mvmult -c mvmult.f90

erf: erf.so

erf.so: erf.f90 constants.mod
	${F2PY} ${F2PYFLAGS} ${INC} ${MKL} -m erf -c erf.f90

constants.mod: Constants.f90
	${F90} ${F90FLAGS} ${INC} ${MKL} -c Constants.f90

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
