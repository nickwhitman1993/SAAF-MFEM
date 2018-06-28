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
   #-DCMAKE_CXX_COMPILER=${CXX} \
   #-DUSE_MPI=OFF \
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
echo "***** Make Broom"
echo "**************************************************************************"
echo ""
make ${PARALLEL_BUILD_FLAG} install
