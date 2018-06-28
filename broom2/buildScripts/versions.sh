#! /bin/sh

# Set your favorite compiler
PARALLEL_BUILD_FLAG="-j16"

# Assuming you're running this from the top of the Broom directory.
SRC_DIR=`pwd`/..

# We do an out-of-source build.
BUILD_DIR=${SRC_DIR}/../buildBroom${BUILDNAME}

SUITESPARSE_VERSION=4.4.4
MFEM_VERSION=3.3-dev
# MFEM_VERSION=3.3
# MFEM_VERSION=3.2-dev
# MFEM_VERSION=3.2
# MFEM_VERSION=cd5e25ebb8
# MFEM_VERSION=5e13c23c94

SUITESPARSE_ROOT=${BUILD_DIR}/TPL
MFEM_ROOT=${BUILD_DIR}/TPL

