// Copyright (c) 2026 Genesis Molecular AI
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.InfoTheory.BitClusterer`, a Python module:
//! clustering bits by their correlation.

use crate::data_structs::ExplicitBitVect;
use crate::sim_div_filters::{ClusterMethod, HierarchicalClusterPicker};
use crate::{Error, Result};

/// Clusters a set of bits by their correlation, and maps fingerprints onto
/// the clusters. Mirrors `BitClusterer.BitClusterer`.
#[derive(Clone, Debug, PartialEq)]
pub struct BitClusterer {
    clusters: Vec<Vec<u32>>,
    bid_list: Vec<u32>,
    n_clusters: usize,
    type_: ClusterMethod,
}

impl BitClusterer {
    /// A clusterer of the bits `id_list` into `n_cluster` clusters, by
    /// Ward's method. Mirrors `BitClusterer(idList, nCluster)`.
    pub fn new(id_list: Vec<u32>, n_cluster: usize) -> BitClusterer {
        BitClusterer::new_with_type_(id_list, n_cluster, ClusterMethod::Ward)
    }

    /// As [`new`](Self::new), by the clustering method `type_`. Mirrors
    /// `BitClusterer(idList, nCluster, type)`.
    pub fn new_with_type_(
        id_list: Vec<u32>,
        n_cluster: usize,
        type_: ClusterMethod,
    ) -> BitClusterer {
        BitClusterer {
            clusters: Vec::new(),
            bid_list: id_list,
            n_clusters: n_cluster,
            type_,
        }
    }

    /// Clusters the bits by a correlation matrix, such as
    /// `BitCorrMatGenerator::get_corr_matrix` returns (the lower triangle,
    /// row by row), taking each distance as the reciprocal of the
    /// correlation. Mirrors `BitClusterer.ClusterBits`.
    ///
    /// Errors as `HierarchicalClusterPicker::cluster`, which refuses the
    /// infinite distance of bits never set together.
    pub fn cluster_bits(&mut self, corr_mat: &[f64]) -> Result<()> {
        let dist_mat: Vec<f64> = corr_mat.iter().map(|c| 1.0 / c).collect();
        let too_many = || Error::Value {
            message: "too many bits or clusters".to_owned(),
        };
        let pkr = HierarchicalClusterPicker::new(self.type_);
        let cls = pkr.cluster(
            &dist_mat,
            u32::try_from(self.bid_list.len()).map_err(|_| too_many())?,
            u32::try_from(self.n_clusters).map_err(|_| too_many())?,
        )?;
        self.clusters = cls
            .into_iter()
            .map(|cl| cl.into_iter().map(|i| self.bid_list[i as usize]).collect())
            .collect();
        Ok(())
    }

    /// Sets the clusters, lists of bits. Mirrors `BitClusterer.SetClusters`.
    ///
    /// `Error::Value` for a number of clusters other than the clusterer's
    /// (Python's assertion).
    pub fn set_clusters(&mut self, clusters: Vec<Vec<u32>>) -> Result<()> {
        if clusters.len() != self.n_clusters {
            return Err(Error::Value {
                message: format!("expected {} clusters", self.n_clusters),
            });
        }
        self.clusters = clusters;
        Ok(())
    }

    /// The clusters, lists of bits. Mirrors `BitClusterer.GetClusters`.
    pub fn get_clusters(&self) -> &[Vec<u32>] {
        &self.clusters
    }

    /// For each cluster, how many of its bits the fingerprint sets. Mirrors
    /// `BitClusterer.MapToClusterScores`.
    ///
    /// `Error::Index` for a bit past the fingerprint.
    pub fn map_to_cluster_scores(&self, fp: &ExplicitBitVect) -> Result<Vec<u32>> {
        let mut scores = vec![0; self.n_clusters];
        for (i, cls) in self.clusters.iter().enumerate() {
            for &bid in cls {
                if fp.get_bit(bid)? {
                    scores[i] += 1;
                }
            }
        }
        Ok(scores)
    }

    /// A fingerprint with a bit for each cluster, set when the fingerprint
    /// sets any of the cluster's bits. Mirrors `BitClusterer.MapToClusterFP`.
    ///
    /// `Error::Index` for a bit past the fingerprint.
    pub fn map_to_cluster_fp(&self, fp: &ExplicitBitVect) -> Result<ExplicitBitVect> {
        let n = u32::try_from(self.n_clusters).map_err(|_| Error::Value {
            message: "too many clusters".to_owned(),
        })?;
        let mut ebv = ExplicitBitVect::new(n);
        for (i, cls) in self.clusters.iter().enumerate() {
            for &bid in cls {
                if fp.get_bit(bid)? {
                    ebv.set_bit(i as u32)?;
                    break;
                }
            }
        }
        Ok(ebv)
    }
}
