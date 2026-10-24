# Copyright (c) 2026 Steven Kearnes
# All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

# Sets <out> to TRUE when <version> belongs to release series <series>
# ("2026.09.1" belongs to "2026.09"), else FALSE.
function(rdkit_cxx_version_supported version series out)
  string(FIND "${version}" "${series}." position)
  if(position EQUAL 0)
    set(${out} TRUE PARENT_SCOPE)
  else()
    set(${out} FALSE PARENT_SCOPE)
  endif()
endfunction()

# Stops configuration unless <version> belongs to release series <series>.
function(rdkit_cxx_check_version version series)
  rdkit_cxx_version_supported("${version}" "${series}" ok)
  if(NOT ok)
    message(FATAL_ERROR
      "This rdkit-cxx crate supports RDKit ${series}.*, but found RDKit "
      "${version} at ${rdkit_DIR}. Install librdkit-dev ${series}.* or point "
      "RDKIT_PREFIX at a matching install.")
  endif()
endfunction()
