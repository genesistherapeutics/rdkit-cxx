# Copyright (c) 2026 Steven Kearnes
# All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

include(${CMAKE_CURRENT_LIST_DIR}/../../wrapper/cmake/RDKitCxxVersion.cmake)
foreach(case "2026.09.1:TRUE" "2026.09.2:TRUE" "2026.09.1b1:TRUE"
             "2026.03.4:FALSE" "2027.03.1:FALSE" "2026.9.1:FALSE")
  string(REPLACE ":" ";" case "${case}")
  list(GET case 0 version)
  list(GET case 1 expected)
  rdkit_cxx_version_supported("${version}" "2026.09" ok)
  if(NOT ok STREQUAL expected)
    message(FATAL_ERROR "${version}: expected ${expected}, got ${ok}")
  endif()
endforeach()
message(STATUS "version check tests passed")
