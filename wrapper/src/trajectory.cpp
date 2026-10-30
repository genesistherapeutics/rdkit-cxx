// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/trajectory.h"

#include <limits>
#include <string>

#include <GraphMol/Trajectory/Snapshot.h>
#include <GraphMol/Trajectory/Trajectory.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

// Both readers size each snapshot as numPoints * 3 in unsigned int, which
// wraps for more points, and then write numPoints * 3 coordinates into the
// smaller array.
void check_num_points(std::uint32_t num_points) {
  if (num_points > std::numeric_limits<unsigned int>::max() / 3) {
    throw ValueErrorException("numPoints * 3 overflows unsigned int");
  }
}

// Reads into a three-dimensional trajectory of num_points points, as Python
// reads into the one it is given, and copies the snapshots' coordinates out,
// num_points * 3 per snapshot.
std::uint32_t read_trajectory(unsigned int (*reader)(const std::string &,
                                                     RDKit::Trajectory &),
                              rust::Slice<const std::uint8_t> path,
                              std::uint32_t num_points,
                              rust::Vec<double> &coords) {
  check_num_points(num_points);
  RDKit::Trajectory trajectory(3, num_points);
  unsigned int count = reader(to_path(path), trajectory);
  coords.clear();
  coords.reserve(trajectory.size() * num_points * 3);
  for (unsigned int i = 0; i < count; ++i) {
    const auto &snapshot = trajectory.getSnapshot(i);
    for (unsigned int p = 0; p < num_points; ++p) {
      auto point = snapshot.getPoint3D(p);
      coords.push_back(point.x);
      coords.push_back(point.y);
      coords.push_back(point.z);
    }
  }
  return count;
}

}  // namespace

std::uint32_t read_amber_trajectory(rust::Slice<const std::uint8_t> path,
                                    std::uint32_t num_points,
                                    rust::Vec<double> &coords) {
  // With no points a snapshot reads nothing, so readAmberTrajectory adds
  // empty snapshots without reading on, forever, once its title line ends in
  // a newline.
  if (num_points == 0) {
    throw ValueErrorException("the trajectory has no points");
  }
  return read_trajectory(RDKit::readAmberTrajectory, path, num_points, coords);
}

std::uint32_t read_gromos_trajectory(rust::Slice<const std::uint8_t> path,
                                     std::uint32_t num_points,
                                     rust::Vec<double> &coords) {
  return read_trajectory(RDKit::readGromosTrajectory, path, num_points, coords);
}

}  // namespace rdkit_cxx
