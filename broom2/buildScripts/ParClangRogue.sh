#! /bin/sh

#use gcc-4.9.2p
#use clang-omp-3.5.0

# Set your favorite compiler
#export MKLROOT=/usr/local/tools/mkl-11.2.0
#export MKLLIB=-Wl,--start-group\ ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a\ ${MKLROOT}/lib/intel64/libmkl_core.a\ ${MKLROOT}/lib/intel64/libmkl_sequential.a\ -Wl,--end-group\ -lpthread\ -ldl
#export CC=gcc-4.6.1
#export MPICC=mpigcc-4.6.1
#export CXX=g++-4.6.1
#export MPICXX=mpig++-4.6.1
export SRC_DIR=$PWD/..
export BLASROOT=${SRC_DIR}/TPL/OpenBLAS
#export BLASLIB=${BLASROOT}/libopenblas.a\ -lpthread\ -ldl
export BLASLIB=-lopenblas
export CC=gcc
export MPICC=mpicc
export CXX=g++
export MPICXX=mpic++
export F77=gfortran
export F90=gfortran
export CFLAGS=-g\ -Wall\ -O3
export EXTRALIBS=-lpthread\ -lm\ -lrt\ -ldl
export CXXFLAGS=${CFLAGS} #\ -fsanitize=undefined,address,integer
export FFLAGS=${CFLAGS}
export F90FLAGS=${FFLAGS}
export OPENMP_FLAGS=''
export BUILDNAME='ParClang'

