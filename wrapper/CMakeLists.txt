# Copyright (c) 2026 Steven Kearnes
# All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

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
  src/mol_list.cpp
  src/mol_bundle.cpp
  src/molops.cpp
  src/substruct.cpp
  src/periodic_table.cpp
  src/monomer_info.cpp
  src/stereo.cpp
  src/substance_group.cpp
  src/query.cpp
  src/resonance.cpp
  src/data_structs.cpp
  src/fpb.cpp
  src/grid.cpp
  src/pickle_checks.cpp
  src/perception.cpp
  src/chirality.cpp
  src/fingerprints.cpp
  src/fingerprint_generator.cpp
  src/adjust_query.cpp
  src/mol_checks.cpp
  src/depiction_lock.cpp
  src/partial_charges.cpp
  src/cip_labeler.cpp
  src/mol_hash.cpp
  src/mhfp.cpp
  src/reduced_graphs.cpp
  src/sln_parse.cpp
  src/conformer_parser.cpp
  src/mmpa.cpp
  src/deprotect.cpp
  src/graph.cpp
  src/transforms.cpp
  src/rdqueries.cpp
  src/descriptors.cpp
  src/descriptors_3d.cpp
  src/dclv.cpp
  src/properties.cpp
  src/standardize.cpp
  src/standardize_fragment.cpp
  src/standardize_charge.cpp
  src/standardize_normalize.cpp
  src/file_checks.cpp
  src/streams.cpp
  src/png.cpp
  src/suppliers.cpp
  src/writers.cpp
)

# RDKit targets the wrappers link, without the RDKit:: prefix (B4).
set(RDKIT_CXX_RDKIT_TARGETS
  RDGeneral
  GraphMol
  SmilesParse
  FileParsers
  MarvinParser
  Depictor
  SubstructMatch
  DataStructs
  RDGeometryLib
  Fingerprints
  GenericGroups
  Subgraphs
  ChemTransforms
  PartialCharges
  CIPLabeler
  MolHash
  ReducedGraphs
  SLNParse
  MMPA
  ChemReactions
  Deprotect
  Descriptors
  MolStandardize
)

add_library(rdkitcxx STATIC ${RDKIT_CXX_SOURCES} ${RDKIT_CXX_GENERATED_SOURCES})
target_include_directories(rdkitcxx PRIVATE
  include
  ${RDKIT_CXX_GENERATED_INCLUDE_DIR}
  ${RDKit_INCLUDE_DIRS}
  # Declarations RDKit's installed headers need (compat/maeparser).
  compat
)
foreach(target IN LISTS RDKIT_CXX_RDKIT_TARGETS)
  target_link_libraries(rdkitcxx PUBLIC RDKit::${target})
  string(APPEND link_libraries "$<TARGET_LINKER_FILE:RDKit::${target}>\n")
endforeach()

# MolBundle's pickles instantiate Boost.Serialization archives in the
# wrappers, when RDKit was built with them.
if(TARGET Boost::serialization)
  target_link_libraries(rdkitcxx PUBLIC Boost::serialization)
  string(APPEND link_libraries "$<TARGET_LINKER_FILE:Boost::serialization>\n")
endif()

# The build script links these, since Cargo, not CMake, links the crate.
file(GENERATE OUTPUT ${CMAKE_BINARY_DIR}/link-libraries.txt
  CONTENT "${link_libraries}")
