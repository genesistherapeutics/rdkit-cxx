# Copyright (c) 2026, Steven Kearnes
# SPDX-License-Identifier: BSD-3-Clause
#
# Builds the C++ half of the rdkit-cxx crate: the hand-written wrapper
# functions and the C++ that cxx generates for each bridge. The crate's build
# script passes:
#   RDKIT_CXX_SUPPORTED_SERIES        RDKit release series, such as 2026.09
#   RDKIT_CXX_GENERATED_SOURCES       ;-separated generated .rs.cc files
#   RDKIT_CXX_GENERATED_INCLUDE_DIR   directory holding rust/cxx.h and the
#                                     generated bridge headers

cmake_minimum_required(VERSION 3.22)
project(rdkitcxx LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(cmake/RDKitCxxVersion.cmake)

find_package(rdkit CONFIG REQUIRED)
rdkit_cxx_check_version("${rdkit_VERSION}" "${RDKIT_CXX_SUPPORTED_SERIES}")

# One wrapper source per area (CONVENTIONS.md B4).
set(RDKIT_CXX_SOURCES
  src/rd_base.cpp
  src/errors.cpp
  src/mol.cpp
  src/molfiles.cpp
  src/atom.cpp
  src/bond.cpp
  src/props.cpp
  src/geometry.cpp
  src/conformer.cpp
  src/ring_info.cpp
)

# RDKit targets the wrappers link, without the RDKit:: prefix (B4).
set(RDKIT_CXX_RDKIT_TARGETS
  RDGeneral
  GraphMol
  SmilesParse
  FileParsers
)

add_library(rdkitcxx STATIC ${RDKIT_CXX_SOURCES} ${RDKIT_CXX_GENERATED_SOURCES})
target_include_directories(rdkitcxx PRIVATE
  include
  ${RDKIT_CXX_GENERATED_INCLUDE_DIR}
  ${RDKit_INCLUDE_DIRS}
)
foreach(target IN LISTS RDKIT_CXX_RDKIT_TARGETS)
  target_link_libraries(rdkitcxx PUBLIC RDKit::${target})
  string(APPEND link_libraries "$<TARGET_LINKER_FILE:RDKit::${target}>\n")
endforeach()

# The build script links these, since Cargo, not CMake, links the crate.
file(GENERATE OUTPUT ${CMAKE_BINARY_DIR}/link-libraries.txt
  CONTENT "${link_libraries}")
