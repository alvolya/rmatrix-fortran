# Rmatrix — build
#   make            build ./Rmatrix
#   make test       run tests/run_tests.py on the build
#   make debug      clean rebuild with bounds checking (-O1 -g -fcheck=bounds)
#   make clean
# GSL is located with pkg-config; override with GSL_PREFIX=/path (headers in
# $(GSL_PREFIX)/include, libraries in $(GSL_PREFIX)/lib).  On Linux the
# LAPACK/BLAS link line is LAPACK_LIBS (default -llapack -lblas, e.g.
# LAPACK_LIBS=-lopenblas); macOS uses the Accelerate framework.

OPTIM         = -O2
COMP          = gfortran
LINKER        = gfortran
F90FLAGS      = -fno-automatic $(OPTIM)
CERNLIB_FLAGS = -I. -DCERNLIB_DOUBLE -DCERNLIB_UNIX -w

# ---------- GSL ----------
GSL_PREFIX ?= /usr/local
ifeq ($(origin GSL_PREFIX),command line)
    GSL_CFLAGS = -I$(GSL_PREFIX)/include
    GSL_LIBS   = -L$(GSL_PREFIX)/lib -lgsl -lgslcblas
else
    GSL_CFLAGS ?= $(shell pkg-config --cflags gsl 2>/dev/null)
    GSL_LIBS   ?= $(shell pkg-config --libs gsl 2>/dev/null)
    ifeq ($(strip $(GSL_LIBS)),)
        GSL_CFLAGS = -I$(GSL_PREFIX)/include
        GSL_LIBS   = -L$(GSL_PREFIX)/lib -lgsl -lgslcblas
    endif
endif

# ---------- platform ----------
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    SDKROOT  := $(shell xcrun --show-sdk-path 2>/dev/null)
    MACOSV   := -mmacosx-version-min=$(shell xcrun --show-sdk-version 2>/dev/null)
    F90FLAGS      += $(MACOSV)
    CERNLIB_FLAGS += $(MACOSV)
    LIBFLS = $(GSL_LIBS) -lm -isysroot $(SDKROOT) -framework Accelerate $(MACOSV)
else
    LAPACK_LIBS ?= -llapack -lblas
    LIBFLS = $(GSL_LIBS) $(LAPACK_LIBS) -lm
endif

# Fortran 90 modules (compiled first, in dependency order)
OBMOD = io_filenames.o channel_data.o resonance_data.o tmat_data.o prep_data.o

# R-matrix physics objects (use modules)
OBRM  = rm.o rmpread.o rmreread.o rmdocalc.o rmblatt.o intgsl.o pretjl.o \
        rmnrphase_complex.o rmcoulph.o rmbbb.o rdef.o apsrb_phase.o \
        rmkinem.o rmbcond.o

# CERNLIB routines compiled from source
OBCERN = wclbes.o \
         r1dp.o r2dp.o r3dp.o r4dp.o r5dp.o r6dp.o r7dp.o r8dp.o \
         clogam64.o cdigam64.o \
         mtlprt.o

.PHONY: all test debug clean

all: Rmatrix

Rmatrix: $(OBMOD) $(OBRM) $(OBCERN)
	$(LINKER) $(OPTIM) -o $@ $(OBMOD) $(OBRM) $(OBCERN) $(LIBFLS)

test: Rmatrix
	python3 tests/run_tests.py --exe ./Rmatrix

debug:
	$(MAKE) clean
	$(MAKE) OPTIM="-O1 -g -fcheck=bounds"

# ---------- Fortran 90 modules ----------
io_filenames.o: io_filenames.f90
	$(COMP) -c $(F90FLAGS) -o $@ $<

channel_data.o: channel_data.f90
	$(COMP) -c $(F90FLAGS) -o $@ $<

resonance_data.o: resonance_data.f90 channel_data.o
	$(COMP) -c $(F90FLAGS) -o $@ $<

tmat_data.o: tmat_data.f90 channel_data.o
	$(COMP) -c $(F90FLAGS) -o $@ $<

prep_data.o: prep_data.f90 channel_data.o
	$(COMP) -c $(F90FLAGS) -o $@ $<

# ---------- Physics routines ----------
rm.o: rm.f90 io_filenames.o
	$(COMP) -c $(F90FLAGS) -o $@ $<

rmpread.o: rmpread.f90 channel_data.o io_filenames.o
	$(COMP) -c $(F90FLAGS) -o $@ $<

rmreread.o: rmreread.f90 channel_data.o resonance_data.o io_filenames.o
	$(COMP) -c $(F90FLAGS) -o $@ $<

rmdocalc.o: rmdocalc.f90 channel_data.o resonance_data.o tmat_data.o prep_data.o io_filenames.o rmkinem.o rmbcond.o
	$(COMP) -c $(F90FLAGS) -o $@ $<

rmblatt.o: rmblatt.f90 channel_data.o resonance_data.o tmat_data.o prep_data.o
	$(COMP) -c $(F90FLAGS) -o $@ $<

intgsl.o: intgsl.f90
	$(COMP) -c $(F90FLAGS) $(GSL_CFLAGS) -o $@ $<

pretjl.o: pretjl.f90 channel_data.o resonance_data.o tmat_data.o prep_data.o intgsl.o
	$(COMP) -c $(F90FLAGS) -o $@ $<

rmnrphase_complex.o: rmnrphase_complex.f90 channel_data.o prep_data.o
	$(COMP) -c $(F90FLAGS) -o $@ $<

rmcoulph.o: rmcoulph.f90 channel_data.o prep_data.o
	$(COMP) -c $(F90FLAGS) -o $@ $<

rmbbb.o: rmbbb.f90 channel_data.o tmat_data.o
	$(COMP) -c $(F90FLAGS) -o $@ $<

rdef.o: rdef.f90 channel_data.o resonance_data.o prep_data.o intgsl.o rmkinem.o
	$(COMP) -c $(F90FLAGS) -o $@ $<

apsrb_phase.o: apsrb_phase.f90 channel_data.o
	$(COMP) -c $(F90FLAGS) -o $@ $<

rmkinem.o: rmkinem.f90 channel_data.o prep_data.o
	$(COMP) -c $(F90FLAGS) -o $@ $<

rmbcond.o: rmbcond.f90 channel_data.o resonance_data.o prep_data.o intgsl.o
	$(COMP) -c $(F90FLAGS) -o $@ $<

# ---------- CERNLIB routines ----------
%.o: %.F
	$(COMP) -c -fno-automatic $(OPTIM) $(CERNLIB_FLAGS) -o $@ $<

mtlprt.o: mtlprt.for
	$(COMP) -c -fno-automatic $(OPTIM) $(CERNLIB_FLAGS) -o $@ $<

clean:
	rm -f *.o *.mod Rmatrix
