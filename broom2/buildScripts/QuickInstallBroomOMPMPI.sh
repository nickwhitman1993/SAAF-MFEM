#! /bin/sh

echo ""
echo "**************************************************************************"
echo "***** Using externally set variables for compiling:"
echo "*****      CC=${CC}"
echo "*****      CXX=${CXX}"
echo "*****      CXXFLAGS=${CXXFLAGS}"
echo "*****      OPENMP_FLAGS=${OPENMP_FLAGS}"
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
echo "***** Compress the source directory. This needs to be done when making changes"
echo "***** to MFEM locally, i.e. mfem-3.2-dev."
echo  "**************************************************************************"
echo ""
echo "DON'T FORGET - we compress the mfem source so it needs to exist..."
cd ${SRC_DIR}/TPL/
tar -jcvf mfem-${MFEM_VERSION}.tar.bz2 mfem-${MFEM_VERSION}
echo "**************************************************************************"
echo ""
cd ${BUILD_DIR}
mkdir -p ${BUILD_DIR}/Broom

echo ""
echo "**************************************************************************"
echo "***** Build the make files.  If you change things like the compilers or options,"
echo "***** cmake requires that you blow away the directory"
echo "**************************************************************************"
echo ""
cd ${BUILD_DIR}/Broom
rm -rf *

echo ""
echo "**************************************************************************"
echo "***** Running CMake"
echo "**************************************************************************"
echo ""
cmake  -DCMAKE_C_COMPILER=${CC} \
   -DCMAKE_CXX_COMPILER=${MPICXX} \
   -DCMAKE_CXX_FLAGS="${CXXFLAGS} ${OPENMP_FLAGS}" \
   -DCMAKE_BUILD_TYPE=None \
   -DUSE_OPENMP=ON \
   -DUSE_MPI=ON \
   -DUSE_METIS=OFF \
   -DUSE_SUITESPARSE=ON \
   -DSUITESPARSE_ROOT="${SUITESPARSE_ROOT}" \
   -DMFEM_ROOT="${MFEM_ROOT}" \
   -DEXTRA_LIBS="${MKLLIB}" \
   -DLINK_BLAS=OFF \
   ${SRC_DIR}

echo ""
echo "**************************************************************************"
echo "***** Uncompressing MFEM"
echo "**************************************************************************"
echo ""
cd ${BUILD_DIR}/TPL
pwd
rm -rf mfem-${MFEM_VERSION}
# tar -jxvf ${SRC_DIR}/TPL/mfem-${MFEM_VERSION}.tar.bz2
tar -vxjf ${SRC_DIR}/TPL/mfem-${MFEM_VERSION}.tar.bz2
ls

echo ""
echo "**************************************************************************"
echo "***** Building MFEM"
echo "**************************************************************************"
echo ""
pwd
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
echo "**************************************************************************"
echo "***** Make Broom"
echo "**************************************************************************"
echo ""
# make ${PARALLEL_BUILD_FLAG} install
cd ${BUILD_DIR}/Broom/src/
make
