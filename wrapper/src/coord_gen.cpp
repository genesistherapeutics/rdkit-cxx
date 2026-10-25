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
#include <GraphMol/ROMol.h>
#include <RDGeneral/Exceptions.h>

#include "rdkit-cxx/depiction_checks.h"
#include "rdkit-cxx/depiction_lock.h"
#include "rdkit-cxx/errors.h"
#include "rdkit-cxx/ffi/coord_gen.rs.h"
#include "rdkit-cxx/util.h"

// GraphMol/CoordGen.h includes coordgen's own headers, which conda-forge's
// librdkit-dev does not install, so this file declares the two pieces of it
// that it uses: CoordGenParams, copied token for token (CMakeLists.txt fails
// the build if the installed header's differs), and addCoords<ROMol>, which
// libRDKitDepictor defines and exports (RDDepictor.cpp instantiates it). The
// copy names coordgen's precision constants, which sketcherMinimizer.h
// defines like these, so both definitions of CoordGenParams are the same
// (the constants are const, internal, and never odr-used).
static const float SKETCHER_STANDARD_PRECISION = 1.f;
static const float SKETCHER_QUICK_PRECISION = 0.2f;
static const float SKETCHER_BEST_PRECISION = 3.f;

namespace RDKit {
namespace CoordGen {

// clang-format off
// BEGIN CoordGen.h
struct CoordGenParams {
  const float sketcherCoarsePrecision = 0.01f;
  const float sketcherStandardPrecision = SKETCHER_STANDARD_PRECISION;
  const float sketcherBestPrecision = SKETCHER_BEST_PRECISION;
  const float sketcherQuickPrecision = SKETCHER_QUICK_PRECISION;
  RDGeom::INT_POINT2D_MAP
      coordMap;  // coordinates for fixing particular atoms of a template
  const ROMol *templateMol = nullptr;  // a molecule to use as a template
  double coordgenScaling = 50.0;       // at the time this was written, coordgen
  // returned coordinates with a single bond
  // length of 50.
  std::string templateFileDir = "";
  float minimizerPrecision =
      sketcherCoarsePrecision;     // controls sketch precision
  bool dbg_useConstrained = true;  // debugging
  bool dbg_useFixed = false;       // debugging
  bool treatNonterminalBondsToMetalAsZeroOrder =
      false;  // set non-terminal bonds to metal atoms to be zero-order bonds
};
// END CoordGen.h
// clang-format on

template <typename T>
unsigned int addCoords(T &mol, const CoordGenParams *params = nullptr);

}  // namespace CoordGen
}  // namespace RDKit

namespace rdkit_cxx {

namespace {

using RDKit::CoordGen::CoordGenParams;

// The parameters add_coords_default uses. RDKit's are a header static,
// CoordGen::defaultParams, one copy per translation unit, so this file keeps
// its own, as rdCoordGen's wrapper does. Read and written under
// lock_depiction().
CoordGenParams &default_params() {
  static CoordGenParams params;
  return params;
}

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

void add_coords_default(RWMol &mol) {
  auto lock = lock_depiction();
  add_coords_with(mol, default_params());
}

// As SetDefaultTemplateFileDir in External/CoordGen/nbWrap/rdCoordGen.cpp.
void set_default_template_file_dir(rust::Slice<const std::uint8_t> dir) {
  auto path = to_path(dir);
  auto lock = lock_depiction();
  default_params().templateFileDir = path;
}

}  // namespace rdkit_cxx
