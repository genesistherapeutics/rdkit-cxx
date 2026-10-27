// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Pins the defaults of the parameter structs of the ports of rdkit.ML,
// rdkit.DataStructs, rdkit.SimDivFilters.SimilarityPickers, and
// rdkit.utils to Python's signatures in RDKit 2024.03.6 (CONVENTIONS.md
// A2).

use rdkit_cxx::data_structs::FoldToTargetDensityParams;
use rdkit_cxx::ml::cluster::clusters::{ClusterParams, PrintParams};
use rdkit_cxx::ml::data::data_utils::{CountResultsParams, FilterDataParams, TextToDataParams};
use rdkit_cxx::ml::data::ml_data::MLDataSetParams;
use rdkit_cxx::ml::data::split_data::SplitIndicesParams;
use rdkit_cxx::ml::files::{DataType, ReFileParams, ReadDataFileParams};
use rdkit_cxx::ml::info_theory::bit_rank::{
    RankBitsParams, calc_info_gains_with_n_possible_bit_vals,
};
use rdkit_cxx::ml::slt::risk::CherkasskyRiskBoundParams;
use rdkit_cxx::sim_div_filters::similarity_pickers::PickerParams;
use rdkit_cxx::utils::chemutils::ConfigToNumElectronsParams;

#[test]
fn defaults_are_pythons() {
    // Cluster(metric=0.0, children=None, position=None, index=-1, name=None, data=None)
    assert_eq!(
        ClusterParams::default(),
        ClusterParams {
            metric: 0.0,
            children: vec![],
            position: vec![],
            index: -1,
            name: None,
            data: None
        }
    );
    // Print(level=0, showData=0, offset='\t')
    assert_eq!(
        PrintParams::default(),
        PrintParams {
            level: 0,
            show_data: false,
            offset: "\t".to_owned()
        }
    );
    // FilterData(inData, val, frac, col=-1, indicesToUse=None, indicesOnly=0)
    assert_eq!(
        FilterDataParams::default(),
        FilterDataParams {
            col: -1,
            indices_to_use: None
        }
    );
    // CountResults(inData, col=-1, bounds=None)
    assert_eq!(
        CountResultsParams::default(),
        CountResultsParams {
            col: -1,
            bounds: None
        }
    );
    // TextToData(reader, ignoreCols=[], onlyCols=None)
    assert_eq!(
        TextToDataParams::default(),
        TextToDataParams {
            ignore_cols: vec![],
            only_cols: None
        }
    );
    // MLDataSet(data, nVars=None, nPts=None, nPossibleVals=None, qBounds=None,
    // varNames=None, ptNames=None, nResults=1)
    assert_eq!(
        MLDataSetParams::default(),
        MLDataSetParams {
            n_vars: None,
            n_pts: None,
            n_possible_vals: None,
            q_bounds: None,
            var_names: None,
            pt_names: None,
            n_results: 1
        }
    );
    // SplitIndices(nPts, frac, silent=1, legacy=0, replacement=0)
    assert_eq!(
        SplitIndicesParams::default(),
        SplitIndicesParams {
            legacy: false,
            replacement: false
        }
    );
    // ReFile(fileName, mode='r', comment='#', trailer='\n')
    assert_eq!(
        ReFileParams::default(),
        ReFileParams {
            comment: "#".to_owned(),
            trailer: Some("\n".to_owned())
        }
    );
    // ReadDataFile(fileName, comment='#', depVarCol=0, dataType=float)
    assert_eq!(
        ReadDataFileParams::default(),
        ReadDataFileParams {
            comment: "#".to_owned(),
            dep_var_col: 0,
            data_type: DataType::Float
        }
    );
    // RankBits(bitVects, actVals, nPossibleBitVals=2, metricFunc=CalcInfoGains)
    let rank = RankBitsParams::default();
    assert_eq!(rank.n_possible_bit_vals, 2);
    let vects = vec![vec![0, 1], vec![1, 1]];
    assert_eq!(
        (rank.metric_func)(&vects, &[0, 1], 2, 2).unwrap(),
        calc_info_gains_with_n_possible_bit_vals(&vects, &[0, 1], 2, 2).unwrap()
    );
    // CherkasskyRiskBound(VCDim, nData, nWrong, conf, a1=1.0, a2=2.0)
    assert_eq!(
        CherkasskyRiskBoundParams::default(),
        CherkasskyRiskBoundParams { a1: 1.0, a2: 2.0 }
    );
    // FoldToTargetDensity(fp, density=0.3, minLength=64)
    assert_eq!(
        FoldToTargetDensityParams::default(),
        FoldToTargetDensityParams {
            density: 0.3,
            min_length: 64
        }
    );
    // TopNOverallPicker(numToPick=10, probeFps=None, dataSet=None,
    // simMetric=TanimotoSimilarity), and SpreadPicker's
    let picker = PickerParams::default();
    assert_eq!(picker.num_to_pick, 10);
    assert!(picker.probe_fps.is_empty() && picker.data_set.is_empty());
    let a = rdkit_cxx::data_structs::ExplicitBitVect::new(8);
    let mut b = a.clone();
    b.set_bit(1).unwrap();
    let mut c = b.clone();
    c.set_bit(2).unwrap();
    assert_eq!((picker.sim_metric)(&b, &c).unwrap(), 0.5);
    // ConfigToNumElectrons(config, ignoreFullD=0, ignoreFullF=0)
    assert_eq!(
        ConfigToNumElectronsParams::default(),
        ConfigToNumElectronsParams {
            ignore_full_d: false,
            ignore_full_f: false
        }
    );
}
