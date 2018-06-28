#! /bin/sh


echo ""
echo "**************************************************************************"
echo "***** Using externally set variables for compiling:"
echo "*****      CC=${CC}"
echo "*****      CXX=${CXX}"
echo "*****      CXXFLAGS=${CXXFLAGS}"
echo "*****      OPENMP_FLAGS=${OPENMP_FLAGS}"
echo "*****      EXTRALIBS=${EXTRALIBS}"
echo "*****      SUITESPARSE_ROOT=${SUITESPARSE_ROOT}"
echo "*****      MKLLIB=${MKLLIB}"
echo "*****      SRC_DIR=${SRC_DIR}"
echo "*****      BUILD_DIR=${BUILD_DIR}"
echo "*****  OPENMP_FLAGS are added onto CXXFLAGS for the appropriate parts"
echo "**************************************************************************"
echo ""

echo ""
echo "**************************************************************************"
echo "***** Making build directories"
echo "***** Source directory: ${SRC_DIR}"
echo "***** Build directory: ${BUILD_DIR}"
echo "**************************************************************************"
echo ""
mkdir -p ${BUILD_DIR}
mkdir -p ${BUILD_DIR}/TPL
mkdir -p ${BUILD_DIR}/TPL/lib
mkdir -p ${BUILD_DIR}/TPL/include

echo ""
echo "**************************************************************************"
echo "***** Uncompressing SuiteSparse"
echo "**************************************************************************"
echo ""
cd ${BUILD_DIR}/TPL
rm -rf SuiteSparse
tar -xvjf ${SRC_DIR}/TPL/SuiteSparse-${SUITESPARSE_VERSION}.tar.bz2
cd SuiteSparse

echo ""
echo "**************************************************************************"
echo "***** Building SuiteSparse"
echo "**************************************************************************"
echo ""
make \
     CC="${CC}" \
     CXX="${MPICXX}" \
     F77="${F77}" \
     CF="${CFLAGS} ${OPENMP_FLAGS}" \
     CXXFLAGS="${CXXFLAGS} ${OPENMP_FLAGS}" \
     INSTALL_LIB="${SUITESPARSE_ROOT}/lib" \
     INSTALL_INCLUDE="${SUITESPARSE_ROOT}/include" \
     BLAS="${MKLLIB}" \
     LAPACK="" \
     LIB="${EXTRALIBS}" \
     CHOLMOD_CONFIG="-DNPARTITION" 
make \
     CC="${CC}" \
     CXX="${MPICXX}" \
     F77="${F77}" \
     CF="${CFLAGS} ${OPENMP_FLAGS}" \
     CXXFLAGS="${CXXFLAGS} ${OPENMP_FLAGS}" \
     INSTALL_LIB="${SUITESPARSE_ROOT}/lib" \
     INSTALL_INCLUDE="${SUITESPARSE_ROOT}/include" \
     BLAS="${MKLLIB}" \
     LAPACK="" \
     LIB="${EXTRALIBS}" \
     CHOLMOD_CONFIG="-DNPARTITION" \
     install
cd ..

echo ""
echo "**************************************************************************"
echo "***** Uncompressing MFEM"
echo "**************************************************************************"
echo ""
cd ${BUILD_DIR}/TPL
rm -rf mfem-${MFEM_VERSION}
tar -xvjf ${SRC_DIR}/TPL/mfem-${MFEM_VERSION}.tar.bz2

echo ""
echo "**************************************************************************"
echo "***** Building MFEM"
echo "**************************************************************************"
echo ""
cd mfem-${MFEM_VERSION}
make config \
     CXX=${CXX} \
     CXXFLAGS="${CXXFLAGS} ${OPENMP_FLAGS}" \
     MFEM_USE_MPI=NO \
     MFEM_DEBUG=NO \
     MFEM_THREAD_SAFE=YES \
     MFEM_USE_SUITESPARSE=YES \
     SUITESPARSE_OPT="-I${SUITESPARSE_ROOT}/include"
     SUITESPARSE_LIB="-L${SUITESPARSE_ROOT}/lib -lklu -lbtf -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig ${MKLLIB}"
make status
make info
make ${PARALLEL_BUILD_FLAG}
make install PREFIX=${MFEM_ROOT}
cd ..

