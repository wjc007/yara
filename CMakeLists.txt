# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/apps/cuda_mapper/CMakeLists.txt
#
# CMakeLists.txt file for cuda_mapper.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_core_apps_cuda_mapper)
message (STATUS "Configuring core/apps/cuda_mapper")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES OpenMP CUDA)
find_package (SeqAn REQUIRED)

# Search for Thrust.
#set(CUDATHRUST_FIND_REQUIRED true)
#find_package (CUDAThrust QUIET)

# Disable CUDA acceleration if CUDA was not found.
if (NOT SEQAN_HAS_CUDA)
  add_definitions (-DCUDA_DISABLED=1)
  message (STATUS "  CUDA disabled: building cuda_mapper without CUDA acceleration.")
endif (NOT SEQAN_HAS_CUDA)

# Warn if OpenMP was not found.
if (NOT SEQAN_HAS_OPENMP)
  message (STATUS "  OpenMP not found: building cuda_mapper without multi-threading.")
endif (NOT SEQAN_HAS_OPENMP)

# ----------------------------------------------------------------------------
# CUDA Setup
# ----------------------------------------------------------------------------

# Workaround for Xcode 4.6 / CUDA 5.x which doesn't find CUDA framework (@rpath).
#if (APPLE)
#    set_target_properties(${TARGET_NAME} PROPERTIES LINK_FLAGS "-F /Library/Frameworks/")
#endif (APPLE)

# Use colornvcc to let Xcode parse nvcc's output.
#if (APPLE)
#    find_program(COLOR_NVCC colornvcc PATHS ${CMAKE_SOURCE_DIR}/util/cmake NO_DEFAULT_PATH)
#    set(CUDA_NVCC_EXECUTABLE ${COLOR_NVCC})
#endif (APPLE)

# Build from Kepler architecture upwards.
#list (APPEND CUDA_NVCC_FLAGS -arch sm_35)

# Enforce nvcc to use at most 32 registers per thread.
#list (APPEND CUDA_NVCC_FLAGS --maxrregcount 32)

# Activate verbose flag in ptxas to show register usage.
#list (APPEND CUDA_NVCC_FLAGS -Xptxas --verbose)

# Output PTX.
#list (APPEND CUDA_NVCC_FLAGS --keep --source-in-ptx)

# ----------------------------------------------------------------------------
# Arch Setup
# ----------------------------------------------------------------------------

if (CMAKE_COMPILER_IS_GNUCXX)
  set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -mpopcnt")
endif (CMAKE_COMPILER_IS_GNUCXX)

set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -march=native")

if (CMAKE_COMPILER_IS_GNUCXX AND APPLE)
  set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -mno-avx")
endif (CMAKE_COMPILER_IS_GNUCXX AND APPLE)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")
list (APPEND CUDA_NVCC_FLAGS "${SEQAN_NVCC_FLAGS}")

# Update the list of file names below if you add source files to your application.
add_executable(cuda_indexer indexer.cpp
                            types.h
                            options.h
                            genome.h
                            genome_index.h)

if (SEQAN_HAS_CUDA)
  cuda_add_executable(cuda_mapper mapper.cpp
                                  mapper.h
                                  mapper.cu
                                  mapper.cuh
                                  tags.h
                                  types.h
                                  options.h
                                  reads.h
                                  genome.h
                                  genome_index.h)
else ()
  add_executable(cuda_mapper mapper.cpp
                             mapper.h
                             tags.h
                             types.h
                             options.h
                             reads.h
                             genome.h
                             genome_index.h)
endif ()

# Add dependencies found by find_package (SeqAn).
target_link_libraries (cuda_indexer ${SEQAN_LIBRARIES})
target_link_libraries (cuda_mapper ${SEQAN_LIBRARIES})

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (cuda_indexer cuda_mapper)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install cuda_mapper in ${PREFIX}/bin directory
install (TARGETS cuda_indexer cuda_mapper
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/cuda_mapper for SeqAn release builds.
install (FILES LICENSE README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (cuda_mapper)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable cuda_mapper in CTD structure.
#set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} cuda_mapper CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:cuda_mapper")
  set (CPACK_PACKAGE_NAME "cuda_mapper")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "CUDA Mapper")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Enrico Siragusa <enrico.siragusa@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (cuda_mapper "cuda_mapper")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:cuda_mapper")
