#! /bin/sh

# Set your favorite compiler
export MKLROOT=/usr/local/tools/mkl-11.2.0
export MKLLIB=-Wl,--start-group\ ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a\ ${MKLROOT}/lib/intel64/libmkl_core.a\ ${MKLROOT}/lib/intel64/libmkl_intel_thread.a\ -Wl,--end-group\ -lpthread\ -ldl
#export MKLLIB=-Wl,--start-group\ ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a\ ${MKLROOT}/lib/intel64/libmkl_core.a\ ${MKLROOT}/lib/intel64/libmkl_sequential.a\ -Wl,--end-group\ -lpthread\ -ldl
export CC=icc-15.0.187
export CXX=icpc-15.0.187
export F77=ifort-15.0.187
export F90=ifort-15.0.187
export EXTRALIBS=-lrt\ -lm
export CFLAGS=-g\ -Wall\ -O3\ -I$(MKLROOT)/include
export CXXFLAGS=${CFLAGS}
export FFLAGS=${CFLAGS}
export F90FLAGS=${FFLAGS}
export OPENMP_FLAGS='-openmp'
export BUILDNAME='OmpIntel'

