# There is a lot of stuff in here that is optional.  Over time, this will be
# documented better, or wrapped with options to make it clearer.  In general,
# if you are having trouble compiling, comment something out.

# Don't allow in-source builds.
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
   message(FATAL_ERROR "In-source builds not allowed. Please make a separate build directory and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

# This may be too restrictive
cmake_minimum_required (VERSION 2.8.7)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

project (Broom)

set (Broom_VERSION_MAJOR 1)
set (Broom_VERSION_MINOR 1)

set (CMAKE_EXPORT_COMPILE_COMMANDS 1)
set (BUILD_SHARED_LIBS OFF)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE PATH "default install path" FORCE )
endif()

# we always want to see everything
if(NOT CMAKE_VERBOSE_MAKEFILE)
  set (CMAKE_VERBOSE_MAKEFILE 1)
endif()

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()

option( LINK_BLAS
   "Optionally run with BLAS and LAPACK libraries, if MFEM was build with them"
   ON )
option( USE_SUITESPARSE
   "Try to find UMFPACK (part of SuiteSparse) for a direct solver in addition to iterative solvers."
   ON )
option( USE_METIS
  "Find Metis for linking with a parallel mfem"
   OFF )
option( USE_HYPRE
  "Find HYPRE for linking with a parallel mfem"
   OFF )
option( USE_OPENMP
   "Build OpenMP components of the code"
   OFF )
option( USE_MPI
  "Build with MPI Support.  Needs Metis, Metis, and a parallel MFEM build"
   OFF )

if( USE_MPI )
   set( BROOM_USE_MPI 1 )
endif()

if( USE_OPENMP )
   set( BROOM_USE_OPENMP 1 )
endif()

include( CheckLibraryExists )
include( CheckIncludeFiles )
include( CheckIncludeFileCXX )

# Now set some default include and library link variables
set( BROOM_INCLUDES )
set( BROOM_LIBS )

set( BROOM_INCLUDES ${BROOM_INCLUDES} ${PROJECT_BINARY_DIR} )

message( STATUS "On ${CMAKE_SYSTEM}, checking for a system-specific high-res timer")
if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
   message( STATUS "Found that we're on a Mac, so we're going to try using the timers in mach/mach_time.h")
   check_include_files( mach/mach_time.h BROOM_HAVE_MACH_TIME )
# This is the wrong way to get this.  Not all ppc's are BG/Q.  This is only
# detecting the front end, not the compute nodes.  We need to fix this for
# later.
elseif( ${CMAKE_SYSTEM} MATCHES "ppc" )
   message( STATUS "Found that we're on an IBM BG/Q, so we're going to try using their high rest timers")
   check_include_files( /bgsys/drivers/ppcfloor/hwi/include/bqc/A2_inlines.h BROOM_HAVE_BGQ_TIME )
elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
   message( STATUS "Found that we're on a regular Linux system, so we're going to try using a the real time library")
   check_library_exists( rt clock_gettime "" BROOM_HAVE_RT_TIME )
   if( BROOM_HAVE_RT_TIME )
      set( BROOM_LIBS ${BROOM_LIBS} "-lrt" )
   endif()
endif()

# Find MFEM, our main library
MESSAGE( STATUS "Looking for MFEM's mfem.hpp" )
find_path( MFEM_INCLUDEDIR
  mfem.hpp
  PATHS ${MFEM_ROOT}/include
  ${MFEM_ROOT}
  ${PROJECT_BINARY_DIR}/../TPL/mfem
  ${PROJECT_SOURCE_DIR}/../TPL/mfem
  ${PROJECT_SOURCE_DIR}/../mfem
  DOC "Path to top of MFEM include directory" )

if( NOT MFEM_INCLUDEDIR )
  message( FATAL_ERROR  "\nMFEM_INCLUDEDIR not found.  MFEM must be installed in directory named 'mfem' next to Broom, specified with a root path '-DMFEM_ROOT=yourpath', or explicity with -DMFEM_INCLUDEDIR=yourpath.")
endif()
MESSAGE( STATUS "Looking for MFEM's mfem.hpp - found in ${MFEM_INCLUDEDIR}" )
set( BROOM_INCLUDES ${BROOM_INCLUDES} ${MFEM_INCLUDEDIR} )
include_directories(${BROOM_INCLUDES})

MESSAGE( STATUS "Looking for MFEM's library " )
find_library(MFEM_LIB
  NAMES mfem
  PATHS ${MFEM_INCLUDEDIR}/lib
  ${MFEM_INCLUDEDIR}
  ${MFEM_ROOT}/lib
  ${MFEM_ROOT}
  ${PROJECT_BINARY_DIR}/../TPL/mfem
  ${PROJECT_SOURCE_DIR}/../TPL/mfem
  ${PROJECT_SOURCE_DIR}/../mfem
  DOC "Path to top of MFEM library directory" )
if( NOT MFEM_LIB )
  message( FATAL_ERROR  "\nMFEM_LIB library not found.  Specify a path with '-DMFEM_LIB=/path/to/libmfem.a' or '-DMFEM_ROOT=yourpath'.")
endif()
MESSAGE( STATUS "Looking for library MFEM_LIB - found ${MFEM_LIB}" )
set( BROOM_LIBS ${BROOM_LIBS} ${MFEM_LIB} )

if(USE_SUITESPARSE)
   # This is NOT how we should do this, but it's quick and get's me building.
   if(SUITESPARSE_ROOT)
     set( UMFPACK_PATH ${SUITESPARSE_ROOT} )
   else()
     set( UMFPACK_PATH ~/local )
   endif()

   find_library( KLU_LIB klu PATHS ${UMFPACK_PATH}/lib NO_DEFAULT_PATH)
   find_library( BTF_LIB btf PATHS ${UMFPACK_PATH}/lib NO_DEFAULT_PATH)
   find_library( UMFPACK_LIB umfpack PATHS ${UMFPACK_PATH}/lib NO_DEFAULT_PATH)
   find_library( AMD_LIB amd PATHS ${UMFPACK_PATH}/lib NO_DEFAULT_PATH)
   find_library( CHOLMOD_LIB cholmod PATHS ${UMFPACK_PATH}/lib NO_DEFAULT_PATH)
   find_library( COLAMD_LIB colamd PATHS ${UMFPACK_PATH}/lib NO_DEFAULT_PATH)
   find_library( SUITESPARSE_LIB suitesparseconfig PATHS ${UMFPACK_PATH}/lib NO_DEFAULT_PATH)
   if( KLU_LIB AND BTF_LIB AND UMFPACK_LIB AND AMD_LIB AND CHOLMOD_LIB AND COLAMD_LIB AND SUITESPARSE_LIB )
      message( STATUS "Found UMFPACK_PATH ${UMFPACK_PATH}" )
      message( STATUS "Found UMFPACK_LIB ${UMFPACK_LIB}" )
      message( STATUS "Found KLU_LIB ${KLU_LIB}" )
      message( STATUS "Found BTF_LIB ${BTF_LIB}" )
      message( STATUS "Found UMFPACK_LIB ${UMFPACK_LIB}" )
      message( STATUS "Found AMD_LIB ${AMD_LIB}" )
      message( STATUS "Found CHOLMOD_LIB ${CHOLMOD_LIB}" )
      message( STATUS "Found COLAMD_LIB ${COLAMD_LIB}" )
      message( STATUS "Found SUITESPARSE_LIB ${SUITESPARSE_LIB}" )
      set( BROOM_LIBS ${BROOM_LIBS} ${KLU_LIB} ${BTF_LIB} ${UMFPACK_LIB} ${AMD_LIB} ${CHOLMOD_LIB} ${COLAMD_LIB} ${SUITESPARSE_LIB})
      set( BROOM_INCLUDES ${BROOM_INCLUDES} ${UMFPACK_PATH}/include )
      set( BROOM_HAVE_UMFPACK 1 )
   endif()
endif()

if(USE_HYPRE)
   # This is NOT how we should do this, but it's quick and get's me building.
   if(HYPRE_ROOT)
     set( HYPRE_PATH ${HYPRE_ROOT} )
   else()
     set( HYPRE_PATH ~/local )
   endif()

   find_library( HYPRE_LIB hypre PATHS ${HYPRE_PATH}/lib )
   if( HYPRE_LIB )
     message( STATUS "Found HYPRE_PATH ${HYPRE_PATH}" )
      message( STATUS "Found HYPRE_LIB ${HYPRE_LIB}" )
      set( BROOM_LIBS ${BROOM_LIBS} ${HYPRE_LIB})
      set( BROOM_INCLUDES ${BROOM_INCLUDES} ${HYPRE_PATH}/include )
      set( BROOM_HAVE_HYPRE 1 )
   endif()
endif()

if(USE_METIS)
   # This is NOT how we should do this, but it's quick and get's me building.
   if(METIS_ROOT)
     set( METIS_PATH ${METIS_ROOT} )
   else()
     set( METIS_PATH ~/local )
   endif()

   find_library( METIS_LIB metis PATHS ${METIS_PATH}/lib )
   if( METIS_LIB )
     message( STATUS "Found METIS_PATH ${METIS_PATH}" )
      message( STATUS "Found METIS_LIB ${METIS_LIB}" )
      set( BROOM_LIBS ${BROOM_LIBS} ${METIS_LIB})
      set( BROOM_INCLUDES ${BROOM_INCLUDES} ${METIS_PATH}/include )
      set( BROOM_HAVE_METIS 1 )
   endif()
endif()

#if( LINK_BLAS )
#   # BLAS and LAPACK are used (optionally) by the innards of MFEM or the
#   # direct solver.
#   # TODO: Take the values from the outside, in buildScripts/version.sh, or
#   # pass this information there.
#   find_package (BLAS)
#   if (BLAS_FOUND)
#      set( BROOM_INCLUDES ${BROOM_INCLUDES} ${BLAS_INCLUDE_DIRS} )
#      set( BROOM_LIBS ${BROOM_LIBS} ${BLAS_LIBRARIES} )
#   endif()
#
#   find_package (LAPACK)
#   if (LAPACK_FOUND)
#      set( BROOM_INCLUDES ${BROOM_INCLUDES} ${LAPACK_INCLUDE_DIRS} )
#      set( BROOM_LIBS ${BROOM_LIBS} ${LAPACK_LIBRARIES} )
#   endif()
#endif()

set( BROOM_LIBS ${BROOM_LIBS} ${EXTRA_LIBS} )



message(STATUS "after TPLs BROOM_LIBS = ${BROOM_LIBS}")
message(STATUS "after TPLs BROOM_INCLUDES = ${BROOM_INCLUDES}")

##############################

configure_file (
   "${PROJECT_SOURCE_DIR}/BroomConfig.hh.in"
   "${PROJECT_BINARY_DIR}/BroomConfig.hh"
   )

##############################
# Configure the documentation.

#find_package(Doxygen)
#if(DOXYGEN_FOUND)
#   configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
#   # or do we want to install the documents somewhere else, say ${CMAKE_INSTALL_PREFIX}
#   add_custom_target(doc
#      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
#      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#      COMMENT "Generating API documentation with Doxygen" VERBATIM
#      )
#endif()

##############################
# Enable testing

ENABLE_TESTING()

# This is to let you run across multiple nodes in a batch job at LLNL.
# It problably isn't as general as it needs to be yet.  If you have trouble,
# just set PRERUN_COMMAND to nothing.
find_path(SRUN srun)

if(SRUN)
   set(PRERUN_COMMAND ${SRUN}/srun -n 1 --exclusive --output=screen.J%J.s%s.%N.n%n.t%t.out)
else()
   set(PRERUN_COMMAND )
endif()

# If you have trouble with the srun stuff on your system, uncomment this.
#set(PRERUN_COMMAND )

##############################
# Now include the rest of the project
add_subdirectory ("${PROJECT_SOURCE_DIR}/src")

