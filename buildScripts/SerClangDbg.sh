#! /bin/sh

# Set your favorite compiler
export MKLLIB='-framework Accelerate'
export CC=clang-mp-3.6
export CXX=clang++-mp-3.6
export F77=gfortran-mp-5
export F90=gfortran-mp-5
export CFLAGS='-g -Wall -O3'
export EXTRALIBS=${MKLLIB}' -lm'
export CXXFLAGS=${CFLAGS}'  -fsanitize=undefined,address,integer -fno-sanitize=vptr'
#export CXXFLAGS=${CFLAGS}' -fsanitize=undefined,address,integer'
export FFLAGS=${CFLAGS}
export F90FLAGS=${FFLAGS}
export OPENMP_FLAGS=''
export BUILDNAME='SerClangDbg'

