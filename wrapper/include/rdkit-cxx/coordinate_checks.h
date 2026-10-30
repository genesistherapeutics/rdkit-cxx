/*
 * Copyright (c) 2026 Genesis Molecular AI
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// Checks of conformers, for RDKit's 3D code that trusts their coordinates.

#pragma once

#include <cmath>
#include <cstddef>
#include <string>

#include <GraphMol/Conformer.h>
#include <RDGeneral/Exceptions.h>

namespace rdkit_cxx {

// RDKit's 3D descriptors and principal axes hand matrices built from
// coordinates (covariance, inertia, gyration, and leverage matrices) to
// Eigen's SVD and eigensolvers, which leave their output uninitialized for
// non-finite input, and RDKit reads it. Coordinates within 1e9 keep those
// matrices, and every distance, finite.
inline void check_finite_coordinates(const RDKit::Conformer &conf,
                                     const std::string &what) {
  for (const auto &pos : conf.getPositions()) {
    for (double v : {pos.x, pos.y, pos.z}) {
      if (!(std::abs(v) <= 1e9)) {
        throw ValueErrorException(what +
                                  " needs finite coordinates within 1e9");
      }
    }
  }
}

// GETAWAY and EEM divide by the distances between atoms and solve the
// matrices they make, which two atoms at one position make infinite or NaN.
// Past 1e-6 those entries stay below about 1e6. Call after
// check_finite_coordinates, which keeps the squared distances finite.
inline void check_distinct_positions(const RDKit::Conformer &conf,
                                     const std::string &what) {
  const auto &pos = conf.getPositions();
  for (std::size_t i = 0; i < pos.size(); ++i) {
    for (std::size_t j = i + 1; j < pos.size(); ++j) {
      if ((pos[i] - pos[j]).lengthSq() < 1e-12) {
        throw ValueErrorException(
            what + " needs atoms at distinct positions, but atoms " +
            std::to_string(i) + " and " + std::to_string(j) +
            " are within 1e-6 of each other");
      }
    }
  }
}

}  // namespace rdkit_cxx
