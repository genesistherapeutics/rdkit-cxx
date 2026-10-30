# Copyright (c) 2026 Genesis Molecular AI, Inc.
# All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

# Sets <out> to the RDKit releases of series <series> that <vendor_dir> holds
# headers for: its subdirectories named for a release, such as 2025.09.6 of
# 2025.09. Each holds copies, from that release, of the headers the wrapper
# needs and RDKit does not install (wrapper/vendor/README.md).
function(rdkit_cxx_supported_releases vendor_dir series out)
  file(GLOB entries RELATIVE "${vendor_dir}" "${vendor_dir}/${series}.*")
  set(releases)
  foreach(entry IN LISTS entries)
    if(IS_DIRECTORY "${vendor_dir}/${entry}")
      list(APPEND releases "${entry}")
    endif()
  endforeach()
  list(SORT releases COMPARE NATURAL)
  set(${out} "${releases}" PARENT_SCOPE)
endfunction()

# Sets <out> to TRUE when <version> is one of the list <releases>, else FALSE.
function(rdkit_cxx_release_supported version releases out)
  list(FIND releases "${version}" index)
  if(index EQUAL -1)
    set(${out} FALSE PARENT_SCOPE)
  else()
    set(${out} TRUE PARENT_SCOPE)
  endif()
endfunction()

# Stops configuration unless <version> is a release of series <series> that
# <vendor_dir> holds headers for. Those headers' layouts can change between
# patch releases, and nothing RDKit installs would catch the mismatch, so
# the wrapper builds only against the releases they came from.
function(rdkit_cxx_check_release version series vendor_dir)
  rdkit_cxx_supported_releases("${vendor_dir}" "${series}" releases)
  if(NOT releases)
    message(FATAL_ERROR "${vendor_dir} holds no RDKit ${series} release.")
  endif()
  rdkit_cxx_release_supported("${version}" "${releases}" ok)
  if(NOT ok)
    list(JOIN releases ", " listed)
    list(GET releases -1 latest)
    message(FATAL_ERROR
      "This rdkit-cxx crate supports RDKit ${listed}, but found RDKit "
      "${version} at ${rdkit_DIR}. The wrapper compiles against copies of "
      "RDKit headers that RDKit does not install, whose layouts can change "
      "between releases, so it builds only against the releases they came "
      "from. Install librdkit-dev ${latest}, point RDKIT_PREFIX at a "
      "supported install, or use an rdkit-cxx release that supports RDKit "
      "${version}.")
  endif()
endfunction()
