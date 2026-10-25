/*
 * Copyright (c) 2026 Steven Kearnes
 * All rights reserved.
 * Use of this source code is governed by a BSD-style
 * license that can be found in the LICENSE file.
 */

// rdMolStandardize's validation and standardization pipeline (Pipeline.h).

#pragma once

#include <memory>

#include <GraphMol/MolStandardize/Pipeline.h>

#include "rust/cxx.h"

namespace rdkit_cxx {

// Shared structs, defined by the bridge (src/ffi/standardize_pipeline.rs).
struct PipelineOptionsData;
struct PipelineResultData;

using Pipeline = RDKit::MolStandardize::Pipeline;
using PipelineStatus = RDKit::MolStandardize::PipelineStatus;
using PipelineStage = RDKit::MolStandardize::PipelineStage;

PipelineOptionsData pipeline_options_default();
// RDKit's standardization step parses the options' metal SMARTS and
// dereferences what SmartsToMol returns, null for SMARTS that do not parse,
// and its MetalDisconnector trusts their shape (standardize_normalize.h's
// check_metal_query), so both are checked here.
std::unique_ptr<Pipeline> pipeline_new(const PipelineOptionsData &options);
// The mol block is checked as the mol-block parser reads it (U12).
PipelineResultData pipeline_run(const Pipeline &pipeline, rust::Str mol_data);

}  // namespace rdkit_cxx
