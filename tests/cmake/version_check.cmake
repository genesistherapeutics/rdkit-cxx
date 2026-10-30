# Copyright (c) 2026 Genesis Molecular AI
# All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

include(${CMAKE_CURRENT_LIST_DIR}/../../wrapper/cmake/RDKitCxxVersion.cmake)

# A release is supported only when it is listed, not merely in the series.
foreach(case "2026.09.1:TRUE" "2026.09.2:FALSE" "2026.09.1b1:FALSE"
             "2026.03.4:FALSE" "2027.03.1:FALSE" "2026.9.1:FALSE")
  string(REPLACE ":" ";" case "${case}")
  list(GET case 0 version)
  list(GET case 1 expected)
  rdkit_cxx_release_supported("${version}" "2026.09.1" ok)
  if(NOT ok STREQUAL expected)
    message(FATAL_ERROR "${version}: expected ${expected}, got ${ok}")
  endif()
endforeach()
rdkit_cxx_release_supported("2026.09.2" "2026.09.1;2026.09.2" ok)
if(NOT ok)
  message(FATAL_ERROR "2026.09.2: expected TRUE among two releases, got ${ok}")
endif()

# wrapper/vendor holds headers for the series' releases, and for nothing else
# that looks like one.
set(vendor ${CMAKE_CURRENT_LIST_DIR}/../../wrapper/vendor)
rdkit_cxx_supported_releases("${vendor}" "2026.09" releases)
if(NOT releases)
  message(FATAL_ERROR "wrapper/vendor holds no RDKit 2026.09 release")
endif()
foreach(release IN LISTS releases)
  if(NOT release MATCHES "^2026\\.09\\.[0-9]+$")
    message(FATAL_ERROR "wrapper/vendor/${release} is not an RDKit release")
  endif()
endforeach()
rdkit_cxx_supported_releases("${vendor}" "2027.03" releases)
if(releases)
  message(FATAL_ERROR "wrapper/vendor holds another series: ${releases}")
endif()
message(STATUS "version check tests passed")
