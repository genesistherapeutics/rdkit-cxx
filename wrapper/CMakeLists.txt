# Copyright (c) 2026 Genesis Molecular AI
# All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

# Builds the C++ half of the rdkit-cxx crate: the hand-written wrapper
# functions and the C++ that cxx generates for each bridge. The crate's build
# script passes:
#   RDKIT_CXX_SUPPORTED_SERIES        RDKit release series, such as 2025.03
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
# vendor/<release> holds the headers below that RDKit does not install, so
# only a release with a directory there is supported.
rdkit_cxx_check_release("${rdkit_VERSION}" "${RDKIT_CXX_SUPPORTED_SERIES}"
  ${CMAKE_CURRENT_SOURCE_DIR}/vendor)

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
  src/chemical_features.cpp
  src/mol_chemical_features.cpp
  src/free_sasa.cpp
  src/mif.cpp
  src/shape_align.cpp
  src/metric_matrix_calc.cpp
  src/murtagh_checks.cpp
  src/cluster.cpp
  src/sim_div_pickers.cpp
  src/info_theory.cpp
  src/abbreviations.cpp
  src/enumerate_stereoisomers.cpp
  src/mol_enumerator.cpp
  src/mol_processing.cpp
  src/mol_interchange.cpp
  src/draw.cpp
  src/draw_checks.cpp
  src/draw_drawer.cpp
  src/synthon_space_search.cpp
  src/inchi.cpp
  src/trajectory.cpp
)

# RDKit targets the wrappers link, without the RDKit:: prefix (B4).
set(RDKIT_CXX_RDKIT_TARGETS
  RDGeneral
  GraphMol
  RingDecomposerLib
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
  ChemicalFeatures
  MolChemicalFeatures
  FreeSASALib
  MolInteractionFields
  PubChemShape
  pubchem_align3d
  hc
  SimDivPickers
  InfoTheory
  Abbreviations
  EnumerateStereoisomers
  MolEnumerator
  MolInterchange
  RDChemDrawLib
  SynthonSpaceSearch
  RDInchiLib
)

# DetermineBonds' Hückel option reads YAeHMOP's parameter table, which
# RDKit has when it was built with YAeHMOP.
if(TARGET RDKit::EHTLib)
  list(APPEND RDKIT_CXX_RDKIT_TARGETS EHTLib)
endif()

add_library(rdkitcxx STATIC ${RDKIT_CXX_SOURCES} ${RDKIT_CXX_GENERATED_SOURCES})
# The release the wrapper compiles against, which errors.cpp compares with the
# one loaded. Only that file needs it, so the others' flags stay the same.
set_source_files_properties(src/errors.cpp PROPERTIES
  COMPILE_DEFINITIONS RDKIT_CXX_RDKIT_VERSION="${rdkit_VERSION}")
# conda-forge's macOS build of RDKit 2025.03.6 has no SLN parser; without
# one, the SLN wrappers refuse every call (sln_parse.cpp), and nothing
# throws the parser's exception (errors.cpp).
if(TARGET RDKit::SLNParse)
  list(APPEND RDKIT_CXX_RDKIT_TARGETS SLNParse)
  set_property(SOURCE src/sln_parse.cpp src/errors.cpp APPEND PROPERTY
    COMPILE_DEFINITIONS RDKIT_CXX_HAVE_SLN_PARSE)
endif()
# GraphMol/CoordGen.h includes coordgen's headers, which RDKit's packages
# do not install; third_party/coordgen holds the version RDKit builds in.
target_include_directories(rdkitcxx PRIVATE
  include
  third_party
  ${RDKIT_CXX_GENERATED_INCLUDE_DIR}
  ${RDKit_INCLUDE_DIRS}
  # Declarations RDKit's installed headers need (compat/maeparser).
  compat
  # RDKit's own copies, from the release found (vendor/<release>), of
  # headers its packages do not install. conda-forge's librdkit-dev 2025.03
  # installs none of ML/InfoTheory's headers (RDKit::InfoTheory's); an
  # install that has them finds its own first, earlier on the path.
  # GraphMol/MolDraw2D holds the
  # drawing code's internal headers (DrawText.h, DrawMol.h, and those they
  # include, and MolDraw2DDetails.h), which no package installs
  # (draw_drawer.cpp).
  vendor/${rdkit_VERSION}
)
# RDKit 2025.03's headers test macros that RDConfig.h defines without
# including it: RingInfo.h declares its ring-family members only under
# RDK_USE_URF. Every source reads RDConfig.h first, so it compiles against
# the classes the library was built with.
target_compile_options(rdkitcxx PRIVATE "SHELL:-include RDGeneral/RDConfig.h")
foreach(target IN LISTS RDKIT_CXX_RDKIT_TARGETS)
  target_link_libraries(rdkitcxx PUBLIC RDKit::${target})
  string(APPEND link_libraries "$<TARGET_LINKER_FILE:RDKit::${target}>\n")
endforeach()

# rdMolProcessing reads gzipped files with Boost.Iostreams' gzip filter, as
# RDKit's gzstream does, when RDKit was built with it.
if(TARGET Boost::iostreams)
  target_link_libraries(rdkitcxx PUBLIC Boost::iostreams)
  string(APPEND link_libraries "$<TARGET_LINKER_FILE:Boost::iostreams>\n")
endif()

# RDKit::MolDraw2D cannot be linked like the targets above: conda-forge's
# package names Cairo::Cairo in its link interface without defining it, and
# FreeType's include directory, which librdkit-dev does not install, in its
# include directories, so CMake refuses it. The wrappers take what they need
# by hand: the library; its compile definitions, which say whether it was
# built with Cairo and FreeType (RDConfig.h, which defines them empty, leaves
# out FreeType), so every wrapper source sees the macros the library was
# built with; and Cairo, whose functions MolDraw2DCairo.h calls inline.
get_target_property(rdkit_cxx_moldraw2d_definitions RDKit::MolDraw2D
  INTERFACE_COMPILE_DEFINITIONS)
if(rdkit_cxx_moldraw2d_definitions)
  foreach(definition IN LISTS rdkit_cxx_moldraw2d_definitions)
    string(REGEX REPLACE "^-D" "" definition "${definition}")
    # RDConfig.h defines these empty; matching it avoids a redefinition.
    if(NOT definition MATCHES "=")
      string(APPEND definition "=")
    endif()
    target_compile_definitions(rdkitcxx PRIVATE "${definition}")
  endforeach()
endif()
string(APPEND link_libraries "$<TARGET_LINKER_FILE:RDKit::MolDraw2D>\n")
if(rdkit_cxx_moldraw2d_definitions MATCHES "RDK_BUILD_CAIRO_SUPPORT")
  get_target_property(rdkit_cxx_moldraw2d_location RDKit::MolDraw2D
    IMPORTED_LOCATION_RELEASE)
  if(NOT rdkit_cxx_moldraw2d_location)
    get_target_property(rdkit_cxx_moldraw2d_location RDKit::MolDraw2D
      IMPORTED_LOCATION)
  endif()
  get_filename_component(rdkit_cxx_rdkit_lib_dir
    "${rdkit_cxx_moldraw2d_location}" DIRECTORY)
  get_filename_component(rdkit_cxx_rdkit_prefix
    "${rdkit_cxx_rdkit_lib_dir}" DIRECTORY)
  find_path(RDKIT_CXX_CAIRO_INCLUDE_DIR cairo.h
    HINTS "${rdkit_cxx_rdkit_prefix}/include" PATH_SUFFIXES cairo REQUIRED)
  find_library(RDKIT_CXX_CAIRO_LIBRARY cairo
    HINTS "${rdkit_cxx_rdkit_lib_dir}" REQUIRED)
  target_include_directories(rdkitcxx PRIVATE "${RDKIT_CXX_CAIRO_INCLUDE_DIR}")
  string(APPEND link_libraries "${RDKIT_CXX_CAIRO_LIBRARY}\n")
endif()

# MolBundle's pickles instantiate Boost.Serialization archives in the
# wrappers, when RDKit was built with them.
if(TARGET Boost::serialization)
  target_link_libraries(rdkitcxx PUBLIC Boost::serialization)
  string(APPEND link_libraries "$<TARGET_LINKER_FILE:Boost::serialization>\n")
endif()

# The build script links these, since Cargo, not CMake, links the crate.
file(GENERATE OUTPUT ${CMAKE_BINARY_DIR}/link-libraries.txt
  CONTENT "${link_libraries}")
# What build.rs reruns for when it changes: the install's version file, which
# an update or reinstall of RDKit replaces, and the libraries it links.
file(GENERATE OUTPUT ${CMAKE_BINARY_DIR}/rdkit-files.txt
  CONTENT "${rdkit_DIR}/rdkit-config-version.cmake\n${link_libraries}")
# The compiler CMake compiled the wrapper with, whose C++ standard library
# build.rs links.
file(GENERATE OUTPUT ${CMAKE_BINARY_DIR}/cxx-compiler.txt
  CONTENT "${CMAKE_CXX_COMPILER}\n")
