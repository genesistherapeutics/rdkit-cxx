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
  src/standardize_validate.cpp
  src/standardize_tautomer.cpp
  src/standardize_pipeline.cpp
  src/file_checks.cpp
  src/streams.cpp
  src/png.cpp
  src/suppliers.cpp
  src/writers.cpp
  src/depiction_checks.cpp
  src/depictor.cpp
  src/coord_gen.cpp
  src/threads.cpp
  src/force_field.cpp
  src/force_field_helpers.cpp
  src/control_c_lock.cpp
  src/embedder.cpp
  src/distance_geometry.cpp
  src/mol_transforms.cpp
  src/alignment.cpp
  src/determine_bonds.cpp
  src/shape_helpers.cpp
  src/mol_align.cpp
  src/chem_reactions.cpp
  src/enumerate.cpp
  src/catalog_checks.cpp
  src/mol_catalog.cpp
  src/frag_catalog.cpp
  src/filter_catalog.cpp
  src/tautomer_query.cpp
  src/generalized_substruct.cpp
  src/rascal_mces.cpp
  src/fmcs.cpp
  src/scaffold_network.cpp
  src/substruct_library.cpp
  src/rgroup.cpp
  src/metric_matrix_calc.cpp
  src/murtagh_checks.cpp
  src/cluster.cpp
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
  coordgen
  MolTransforms
  ForceField
  ForceFieldHelpers
  Trajectory
  DistGeometry
  DistGeomHelpers
  Alignment
  DetermineBonds
  ShapeHelpers
  MolAlign
  O3AAlign
  FilterCatalog
  MolStandardize
  Catalogs
  MolCatalog
  FragCatalog
  TautomerQuery
  GeneralizedSubstruct
  RascalMCES
  FMCS
  ScaffoldNetwork
  SubstructLibrary
  RGroupDecomposition
  hc
)

# DetermineBonds' Hückel option reads YAeHMOP's parameter table, which
# RDKit has when it was built with YAeHMOP.
if(TARGET RDKit::EHTLib)
  list(APPEND RDKIT_CXX_RDKIT_TARGETS EHTLib)
endif()

add_library(rdkitcxx STATIC ${RDKIT_CXX_SOURCES} ${RDKIT_CXX_GENERATED_SOURCES})
# GraphMol/CoordGen.h includes coordgen's headers, which RDKit's packages
# do not install; third_party/coordgen holds the version RDKit builds in.
target_include_directories(rdkitcxx PRIVATE
  include
  third_party
  ${RDKIT_CXX_GENERATED_INCLUDE_DIR}
  ${RDKit_INCLUDE_DIRS}
  # Declarations RDKit's installed headers need (compat/maeparser).
  compat
  # conda-forge's librdkit-dev 2026.09 installs DistGeomHelpers/
  # BoundsMatrixBuilder.h without the BoundsMatrixBuilderDetails.h it
  # includes; vendor/ holds RDKit's own copy, from Release_2026_09_1. An
  # install that has it finds its own first, next to the including header.
  vendor
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
