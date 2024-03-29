# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/apps/yara/CMakeLists.txt
#
# CMakeLists.txt file for yara.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.9)
project (seqan_core_apps_yara)
message (STATUS "Configuring core/apps/yara")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search C++11.
find_package(CXX11 REQUIRED)

# Stop here if C++11 was not found.
if (NOT CXX11_FOUND)
    message (STATUS "  C++11 compiler required to build yara.")
    return ()
endif (NOT CXX11_FOUND)

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES OpenMP ZLIB BZIP2) #CUDA
find_package (SeqAn REQUIRED)

# Disable CUDA acceleration if CUDA was not found.
#if (NOT SEQAN_HAS_CUDA)
#  add_definitions (-DCUDA_DISABLED=1)
#  message (STATUS "  CUDA disabled: building yara without CUDA acceleration.")
#endif (NOT SEQAN_HAS_CUDA)

# Warn if OpenMP was not found.
if (NOT SEQAN_HAS_OPENMP)
  message (STATUS "  OpenMP not found: building yara without multi-threading.")
endif (NOT SEQAN_HAS_OPENMP)

# Warn if Zlib was not found.
if (NOT SEQAN_HAS_ZLIB)
  message (STATUS "  Zlib not found: building yara without gzipped input and BAM output.")
endif (NOT SEQAN_HAS_ZLIB)

# Warn if Bzip2 was not found.
if (NOT SEQAN_HAS_BZIP2)
  message (STATUS "  Bzip2 not found: building yara without bzipped input.")
endif (NOT SEQAN_HAS_BZIP2)

# ----------------------------------------------------------------------------
# CUDA Setup
# ----------------------------------------------------------------------------

# Setup CUDA flags.
#seqan_setup_cuda_vars (ARCH sm_35)

# Enforce nvcc to use at most 32 registers per thread.
#list (APPEND CUDA_NVCC_FLAGS --maxrregcount 32)

# Activate verbose flag in ptxas to show register usage.
#list (APPEND CUDA_NVCC_FLAGS -Xptxas --verbose)

# Output PTX.
#list (APPEND CUDA_NVCC_FLAGS --keep --source-in-ptx)

# Workaround for Xcode 4.6 / CUDA 5.x which doesn't find CUDA framework (@rpath).
#if (APPLE)
#    set_target_properties(${TARGET_NAME} PROPERTIES LINK_FLAGS "-F /Library/Frameworks/")
#endif (APPLE)

# ----------------------------------------------------------------------------
# Arch Setup
# ----------------------------------------------------------------------------

if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
  set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -march=native")
endif (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)

if (CMAKE_COMPILER_IS_GNUCXX)
  set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -mpopcnt")
endif (CMAKE_COMPILER_IS_GNUCXX)

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
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${CXX11_CXX_FLAGS}")

# Update the list of file names below if you add source files to your application.
add_executable(yara_indexer indexer.cpp
                            store_genome.h
                            misc_timer.h
                            misc_options.h
                            misc_types.h
                            index_fm.h)

#if (SEQAN_HAS_CUDA)
#  cuda_add_executable(yara_mapper mapper.cpp
#                                  mapper.h
#                                  mapper.cu
#                                  mapper.cuh
#                                  mapper_collector.h
#                                  mapper_classifier.h
#                                  mapper_ranker.h
#                                  mapper_filter.h
#                                  mapper_extender.h
#                                  mapper_verifier.h
#                                  mapper_selector.h
#                                  mapper_aligner.h
#                                  mapper_writer.h
#                                  store_reads.h
#                                  store_genome.h
#                                  misc_tags.h
#                                  misc_timer.h
#                                  misc_options.h
#                                  misc_types.h
#                                  bits_hits.h
#                                  bits_matches.h
#                                  bits_context.h
#                                  bits_seeds.h
#                                  find_extender.h
#                                  find_verifier.h
#                                  index_fm.h)
#else ()
  add_executable(yara_mapper      mapper.cpp
                                  mapper.h
                                  mapper_collector.h
                                  mapper_classifier.h
                                  mapper_ranker.h
                                  mapper_filter.h
                                  mapper_extender.h
                                  mapper_verifier.h
                                  mapper_selector.h
                                  mapper_aligner.h
                                  mapper_writer.h
                                  store_reads.h
                                  store_genome.h
                                  misc_tags.h
                                  misc_timer.h
                                  misc_options.h
                                  misc_types.h
                                  bits_hits.h
                                  bits_matches.h
                                  bits_context.h
                                  bits_seeds.h
                                  find_extender.h
                                  find_verifier.h
                                  index_fm.h)
#endif ()

# Add dependencies found by find_package (SeqAn).
target_link_libraries (yara_indexer ${SEQAN_LIBRARIES})
target_link_libraries (yara_mapper ${SEQAN_LIBRARIES})

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (yara_indexer yara_mapper)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install yara in ${PREFIX}/bin directory
install (TARGETS yara_indexer yara_mapper
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/yara for SeqAn release builds.
install (FILES LICENSE README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (yara_mapper)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable yara_mapper in CTD structure.
#set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} yara_mapper CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:yara")
  set (CPACK_PACKAGE_NAME "yara")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Yara - Yet Another Read Aligner")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Enrico Siragusa <enrico.siragusa@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (yara "yara")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:yara")
