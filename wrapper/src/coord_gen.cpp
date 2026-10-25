// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "rdkit-cxx/coord_gen.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <string>

#include <Geometry/point.h>
// Includes coordgen's headers, from third_party/coordgen.
#include <GraphMol/CoordGen.h>
#include <GraphMol/ROMol.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/depiction_checks.h"
#include "rdkit-cxx/depiction_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/coord_gen.rs.h"
#include "rdkit-cxx/util.h"

namespace rdkit_cxx {

namespace {

using RDKit::CoordGen::CoordGenParams;

// A template coordinate as CoordGen receives it, times coordgenScaling.
void check_template_coordinate(double v, double scaling) {
  if (!(std::abs(v * scaling) <= kMaxCoordGenCoordinate)) {
    throw ValueErrorException(
        "CoordGen template coordinates times coordgenScaling must be finite "
        "and within 5e6");
  }
}

// What CoordGen::addCoords trusts. The caller holds lock_depiction().
void check_params(const RDKit::ROMol &mol, const CoordGenParams &params) {
  // CoordgenMinimizer::growSolutions casts 6 * minimizerPrecision to int,
  // and a negative one rejects every solution, so findBestSolution
  // dereferences the first of none.
  const float precision = params.minimizerPrecision;
  if (!(precision >= 0.f && 6 * precision < 2147483648.f)) {
    throw ValueErrorException("minimizerPrecision must be between 0 and 3.5e8");
  }
  for (const auto &[idx, p] : params.coordMap) {
    check_template_coordinate(p.x, params.coordgenScaling);
    check_template_coordinate(p.y, params.coordgenScaling);
  }
  // addCoords reads a template's coordinates only when it has one
  // conformer.
  if (params.templateMol && params.templateMol->getNumConformers() == 1) {
    for (const auto &p : params.templateMol->getConformer().getPositions()) {
      check_template_coordinate(p.x, params.coordgenScaling);
      check_template_coordinate(p.y, params.coordgenScaling);
    }
  }
  check_coord_gen_mol(mol, params.treatNonterminalBondsToMetalAsZeroOrder);
}

void add_coords_with(RWMol &mol, const CoordGenParams &params) {
  check_params(mol, params);
  RDKit::CoordGen::addCoords<RDKit::ROMol>(mol, &params);
}

}  // namespace

CoordGenParamsData coord_gen_params_default() {
  CoordGenParams params;
  return CoordGenParamsData{
      {},
      {},
      params.coordgenScaling,
      to_bytes(params.templateFileDir),
      params.minimizerPrecision,
      params.dbg_useConstrained,
      params.dbg_useFixed,
      params.treatNonterminalBondsToMetalAsZeroOrder,
      params.sketcherCoarsePrecision,
      params.sketcherStandardPrecision,
      params.sketcherBestPrecision,
      params.sketcherQuickPrecision,
  };
}

// As addCoordsHelper in External/CoordGen/nbWrap/rdCoordGen.cpp. Python's
// SetCoordMap takes any index; CoordGen looks up only the molecule's atoms.
void add_coords(RWMol &mol, const CoordGenParamsData &data,
                const RWMol &template_mol, bool has_template_mol) {
  CoordGenParams params;
  const auto n =
      std::min(data.coord_map_keys.size(), data.coord_map_values.size());
  for (std::size_t i = 0; i < n; ++i) {
    const auto &p = data.coord_map_values[i];
    params.coordMap[static_cast<int>(data.coord_map_keys[i])] =
        RDGeom::Point2D(p.x, p.y);
  }
  params.templateMol = has_template_mol ? &template_mol : nullptr;
  params.coordgenScaling = data.coordgen_scaling;
  params.templateFileDir = to_path(rust::Slice<const std::uint8_t>(
      data.template_file_dir.data(), data.template_file_dir.size()));
  params.minimizerPrecision = data.minimizer_precision;
  params.dbg_useConstrained = data.dbg_use_constrained;
  params.dbg_useFixed = data.dbg_use_fixed;
  params.treatNonterminalBondsToMetalAsZeroOrder =
      data.treat_nonterminal_bonds_to_metal_as_zobs;
  auto lock = lock_depiction();
  add_coords_with(mol, params);
}

// CoordGen::defaultParams is a header static, one copy per translation unit;
// this one is the wrapper's, as rdCoordGen's is its own. Read and written
// under lock_depiction().
void add_coords_default(RWMol &mol) {
  auto lock = lock_depiction();
  add_coords_with(mol, RDKit::CoordGen::defaultParams);
}

// As SetDefaultTemplateFileDir in External/CoordGen/nbWrap/rdCoordGen.cpp.
void set_default_template_file_dir(rust::Slice<const std::uint8_t> dir) {
  auto path = to_path(dir);
  auto lock = lock_depiction();
  RDKit::CoordGen::defaultParams.templateFileDir = path;
}

}  // namespace rdkit_cxx
