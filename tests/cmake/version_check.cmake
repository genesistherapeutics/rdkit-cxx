# Copyright (c) 2026 Steven Kearnes
# All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

include(${CMAKE_CURRENT_LIST_DIR}/../../wrapper/cmake/RDKitCxxVersion.cmake)

# A release is supported only when it is listed, not merely in the series.
foreach(case "2025.09.6:TRUE" "2025.09.7:FALSE" "2025.09.6b1:FALSE"
             "2025.03.4:FALSE" "2026.03.6:FALSE" "2025.9.6:FALSE")
  string(REPLACE ":" ";" case "${case}")
  list(GET case 0 version)
  list(GET case 1 expected)
  rdkit_cxx_release_supported("${version}" "2025.09.6" ok)
  if(NOT ok STREQUAL expected)
    message(FATAL_ERROR "${version}: expected ${expected}, got ${ok}")
  endif()
endforeach()
rdkit_cxx_release_supported("2025.09.7" "2025.09.6;2025.09.7" ok)
if(NOT ok)
  message(FATAL_ERROR "2025.09.7: expected TRUE among two releases, got ${ok}")
endif()

# wrapper/vendor holds headers for the series' releases, and for nothing else
# that looks like one.
set(vendor ${CMAKE_CURRENT_LIST_DIR}/../../wrapper/vendor)
rdkit_cxx_supported_releases("${vendor}" "2025.09" releases)
if(NOT releases)
  message(FATAL_ERROR "wrapper/vendor holds no RDKit 2025.09 release")
endif()
foreach(release IN LISTS releases)
  if(NOT release MATCHES "^2025\\.09\\.[0-9]+$")
    message(FATAL_ERROR "wrapper/vendor/${release} is not an RDKit release")
  endif()
endforeach()
rdkit_cxx_supported_releases("${vendor}" "2026.03" releases)
if(releases)
  message(FATAL_ERROR "wrapper/vendor holds another series: ${releases}")
endif()
message(STATUS "version check tests passed")
