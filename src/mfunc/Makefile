# The compiler(s)
LIFX := $(shell which ifx 2>/dev/null)
ifeq ($(LIFX), )
	F90 = gfortran
else
	F90 = ifx
	F90FLAGS += -heap-arrays
endif

# F90  = gfortran

# Release 
F2PY = FC="${F90}"  python -m numpy.f2py
F90FLAGS += -lpthread 

ifeq ($(MKLROOT), )
	# No MKL
	F90FLAGS += -lblas
else
	F90FLAGS += -lirc -limf -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lmkl_rt -lm -ldl
	INC += -I$(MKLROOT)/include
endif

# NOTE: ifort somehow fails sanity chekc of f2py with meson backend 
F2PYFLAGS =  --backend meson --f90flags=${F90FLAGS}

INC += -I$(CURDIR)
# =============================================================================

# =============================================================================
# Building rules
# =============================================================================

# "make" builds all
all: cleanall erf mvmult

# General rule for building executable from objects. 
# $@ is the name of the target (in this case the executable)
mvmult: mvmult.so

mvmult.so: mvmult.f90 constants.mod
	${F2PY} ${F2PYFLAGS} ${INC}  -m mvmult -c mvmult.f90

erf: erf.so

erf.so: erf.f90 constants.mod
	${F2PY} ${F2PYFLAGS} ${INC} -m erf -c erf.f90

constants.mod: Constants.f90
	${F90} ${F90FLAGS} ${INC} -c Constants.f90

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
