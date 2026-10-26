// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.Pharm2D.SigFactory`, a Python module: the
//! parameters of 2D pharmacophore signatures and the numbering of their
//! bits.

use std::collections::BTreeMap;

use crate::chem::Mol;
use crate::chem::pharm_2d::utils::{
    count_up_to, get_index_combinations, get_possible_scaffolds_with_use_triangle_inequality,
    num_combinations, order_triangle, point_dists,
};
use crate::chem::rd_mol_chemical_features::{GetMolFeatureParams, MolChemicalFeatureFactory};
use crate::data_structs::{IntSparseIntVect, LongSparseIntVect, SparseBitVect};
use crate::{Error, Result};

/// A 2D pharmacophore signature: bits, or counts when the factory uses
/// counts, in a vector with 64-bit indices from 2^31 bits on, as Python
/// chooses its class.
#[derive(Debug)]
pub enum Signature {
    /// Which pharmacophores the molecule has.
    Bits(SparseBitVect),
    /// How many times the molecule has each pharmacophore.
    Counts(IntSparseIntVect),
    /// Counts, for signatures of 2^31 bits or more.
    LongCounts(LongSparseIntVect),
}

/// The molecule's features of `family`, each as its atoms, as Python's
/// `MolChemicalFeatureFactory.GetFeaturesForMol(mol, includeOnly=family)`.
pub(crate) fn feature_atoms(
    factory: &MolChemicalFeatureFactory,
    mol: &Mol,
    family: &str,
) -> Result<Vec<Vec<u32>>> {
    let count = factory.get_num_mol_features_with_include_only(mol, family)?;
    let params = GetMolFeatureParams {
        include_only: family.to_owned(),
        conf_id: -1,
    };
    (0..count)
        .map(|i| {
            Ok(factory
                .get_mol_feature_with_params(mol, i, &params)?
                .get_atom_ids())
        })
        .collect()
}

/// Parameters for [`SigFactory::new_with_params`], with Python's defaults
/// (A2).
#[derive(Clone, Debug, PartialEq)]
pub struct SigFactoryParams {
    /// Counts each pharmacophore instead of setting its bit.
    pub use_counts: bool,
    /// The fewest points in a pharmacophore.
    pub min_point_count: usize,
    /// The most points in a pharmacophore.
    pub max_point_count: usize,
    /// Uses only the shortest path between features; the only kind of
    /// signature Python generates.
    pub shortest_paths_only: bool,
    /// Weights the topological distances by bond order.
    pub include_bond_order: bool,
    /// Feature families to leave out.
    pub skip_feats: Vec<String>,
    /// Leaves out scaffolds whose bins cannot make triangles.
    pub triangle_prune_bins: bool,
}

impl Default for SigFactoryParams {
    fn default() -> SigFactoryParams {
        SigFactoryParams {
            use_counts: false,
            min_point_count: 2,
            max_point_count: 3,
            shortest_paths_only: true,
            include_bond_order: false,
            skip_feats: Vec::new(),
            triangle_prune_bins: true,
        }
    }
}

/// Makes 2D pharmacophore signatures: set the parameters and the distance
/// bins, call [`init`](SigFactory::init), then fingerprint molecules with
/// [`gen_2d_fingerprint`](crate::chem::pharm_2d::generate::gen_2d_fingerprint).
/// As in Python, the parameters are public, and changing one needs another
/// `init`. Mirrors `SigFactory.SigFactory`.
#[derive(Debug)]
pub struct SigFactory {
    /// The factory whose feature families are the pharmacophores' points.
    pub feat_factory: MolChemicalFeatureFactory,
    /// As in [`SigFactoryParams`].
    pub use_counts: bool,
    /// As in [`SigFactoryParams`].
    pub min_point_count: usize,
    /// As in [`SigFactoryParams`].
    pub max_point_count: usize,
    /// As in [`SigFactoryParams`].
    pub shortest_paths_only: bool,
    /// As in [`SigFactoryParams`].
    pub include_bond_order: bool,
    /// As in [`SigFactoryParams`].
    pub skip_feats: Vec<String>,
    /// As in [`SigFactoryParams`].
    pub triangle_prune_bins: bool,
    bins: Vec<(u32, u32)>,
    scaffolds: BTreeMap<usize, Vec<Vec<usize>>>,
    starts: BTreeMap<usize, u64>,
    n_feats: usize,
    sig_size: u64,
}

impl SigFactory {
    /// A factory with Python's defaults and no bins. Mirrors
    /// `SigFactory(featFactory)`.
    pub fn new(feat_factory: MolChemicalFeatureFactory) -> SigFactory {
        SigFactory::new_with_params(feat_factory, &SigFactoryParams::default())
    }

    /// Mirrors `SigFactory(featFactory, useCounts, minPointCount,
    /// maxPointCount, shortestPathsOnly, includeBondOrder, skipFeats,
    /// trianglePruneBins)`.
    pub fn new_with_params(
        feat_factory: MolChemicalFeatureFactory,
        params: &SigFactoryParams,
    ) -> SigFactory {
        SigFactory {
            feat_factory,
            use_counts: params.use_counts,
            min_point_count: params.min_point_count,
            max_point_count: params.max_point_count,
            shortest_paths_only: params.shortest_paths_only,
            include_bond_order: params.include_bond_order,
            skip_feats: params.skip_feats.clone(),
            triangle_prune_bins: params.triangle_prune_bins,
            bins: Vec::new(),
            scaffolds: BTreeMap::new(),
            starts: BTreeMap::new(),
            n_feats: 0,
            sig_size: 0,
        }
    }

    /// Sets the distance bins, each a lower bound and an upper bound past
    /// it, in order, and calls [`init`](SigFactory::init). Mirrors
    /// `SigFactory.SetBins`.
    pub fn set_bins(&mut self, bins: &[(u32, u32)]) -> Result<()> {
        self.bins = bins.to_vec();
        self.init()
    }

    /// Mirrors `SigFactory.GetBins`.
    pub fn get_bins(&self) -> &[(u32, u32)] {
        &self.bins
    }

    /// Mirrors `SigFactory.GetNumBins`.
    pub fn get_num_bins(&self) -> usize {
        self.bins.len()
    }

    /// An empty signature of [`get_sig_size`](SigFactory::get_sig_size)
    /// bits, counts if the factory uses them. Mirrors
    /// `SigFactory.GetSignature`.
    pub fn get_signature(&self) -> Result<Signature> {
        if !self.use_counts {
            let size = u32::try_from(self.sig_size).map_err(|_| Error::Value {
                message: format!("a signature of {} bits is too large", self.sig_size),
            })?;
            Ok(Signature::Bits(SparseBitVect::new(size)))
        } else if self.sig_size < 1 << 31 {
            Ok(Signature::Counts(IntSparseIntVect::new(
                self.sig_size as i32,
            )?))
        } else {
            let size = i64::try_from(self.sig_size).map_err(|_| Error::Value {
                message: format!("a signature of {} bits is too large", self.sig_size),
            })?;
            Ok(Signature::LongCounts(LongSparseIntVect::new(size)?))
        }
    }

    /// The bit's labels and its distances between each pair of points as
    /// bin indices, as Python's `_GetBitSummaryData`.
    fn bit_summary_data(&self, bit_idx: u64) -> Result<(Vec<String>, Vec<Vec<usize>>)> {
        let (n_pts, combo, scaffold) = self.get_bit_info(bit_idx)?;
        let fams = self.get_feat_families();
        let labels = combo.iter().map(|&x| fams[x].clone()).collect();
        let mut d_mat = vec![vec![0; n_pts]; n_pts];
        for (&(i, j), &d) in point_dists(n_pts)?.iter().zip(&scaffold) {
            d_mat[i][j] = d;
            d_mat[j][i] = d;
        }
        Ok((labels, d_mat))
    }

    /// A description of the bit: its features' families, then its matrix of
    /// distance bin indices, row by row (`"Acceptor Donor |0 2|2 0|"`).
    /// Mirrors `SigFactory.GetBitDescription`.
    pub fn get_bit_description(&self, bit_idx: u64) -> Result<String> {
        let (labels, d_mat) = self.bit_summary_data(bit_idx)?;
        let mut res = labels.join(" ") + " ";
        for row in d_mat {
            let row: Vec<String> = row.iter().map(usize::to_string).collect();
            res += &format!("|{}", row.join(" "));
        }
        res += "|";
        Ok(res)
    }

    /// The bin index of each distance, then the scaffold's position among
    /// `scaffolds`, as Python's `_findBinIdx`; `None` for a distance in no
    /// bin or a scaffold not among them.
    fn find_bin_idx(&self, dists: &[u32], scaffolds: &[Vec<usize>]) -> Option<usize> {
        let mut which_bins = Vec::with_capacity(dists.len());
        for &dist in dists {
            let (mut start, mut end) = (0, self.bins.len());
            let mut found = None;
            while start < end {
                let mid = (start + end) / 2;
                let (beg_bin, end_bin) = self.bins[mid];
                if dist < beg_bin {
                    end = mid;
                } else if dist >= end_bin {
                    start = mid + 1;
                } else {
                    found = Some(mid);
                    break;
                }
            }
            which_bins.push(found?);
        }
        scaffolds.iter().position(|s| *s == which_bins)
    }

    /// The families of the factory's features, less those skipped, sorted.
    /// Mirrors `SigFactory.GetFeatFamilies`.
    pub fn get_feat_families(&self) -> Vec<String> {
        let mut fams: Vec<String> = self
            .feat_factory
            .get_feature_families()
            .into_iter()
            .filter(|fam| !self.skip_feats.contains(fam))
            .collect();
        fams.sort();
        fams
    }

    /// The atoms of each of the molecule's features, by family in
    /// [`get_feat_families`](SigFactory::get_feat_families)' order. Mirrors
    /// `SigFactory.GetMolFeats`.
    pub fn get_mol_feats(&self, mol: &Mol) -> Result<Vec<Vec<Vec<u32>>>> {
        self.get_feat_families()
            .iter()
            .map(|fam| feature_atoms(&self.feat_factory, mol, fam))
            .collect()
    }

    /// The bit of the pharmacophore with these features (indices into
    /// [`get_feat_families`](SigFactory::get_feat_families), sorted first)
    /// and distances (in the order of
    /// [`n_point_dist_dict`](crate::chem::pharm_2d::utils::n_point_dist_dict)).
    /// More than 3 points are `Error::Runtime` (Python's
    /// `NotImplementedError`); a number of points outside the factory's,
    /// a feature index past its families, and distances outside its bins
    /// or scaffolds are `Error::Index`; and distances of another number than
    /// the points need are `Error::Value`. Mirrors `SigFactory.GetBitIdx`.
    pub fn get_bit_idx(&self, feat_indices: &[usize], dists: &[u32]) -> Result<u64> {
        self.get_bit_idx_with_sort_indices(feat_indices, dists, true)
    }

    /// Takes the feature indices in order unless `sort_indices`. Mirrors
    /// `SigFactory.GetBitIdx(featIndices, dists, sortIndices)`.
    pub fn get_bit_idx_with_sort_indices(
        &self,
        feat_indices: &[usize],
        dists: &[u32],
        sort_indices: bool,
    ) -> Result<u64> {
        let n_points = feat_indices.len();
        if n_points > 3 {
            return Err(Error::Runtime {
                message: ">3 points not supported".to_owned(),
            });
        }
        if n_points < self.min_point_count || n_points > self.max_point_count {
            return Err(Error::Index {
                index: n_points as u32,
                len: self.max_point_count as u32 + 1,
            });
        }
        let start_idx = *self.starts.get(&n_points).ok_or_else(|| Error::Value {
            message: "the factory needs init".to_owned(),
        })?;
        let mut feat_indices = feat_indices.to_vec();
        if sort_indices {
            feat_indices.sort();
        }
        let max = *feat_indices
            .iter()
            .max()
            .expect("at least min_point_count points");
        if max >= self.n_feats {
            return Err(Error::Index {
                index: max as u32,
                len: self.n_feats as u32,
            });
        }
        let mut dists = dists.to_vec();
        if n_points == 3 {
            (feat_indices, dists) = order_triangle(&feat_indices, &dists)?;
        }
        let scaffolds = self
            .scaffolds
            .get(&dists.len())
            .ok_or_else(|| Error::Value {
                message: format!("{} distances for {n_points} points", dists.len()),
            })?;
        let offset = count_up_to(self.n_feats, n_points, &feat_indices)? * scaffolds.len() as u64;
        // Python raises IndexError for a scaffold its bins lack, and
        // TypeError for a distance in no bin.
        let not_found = || Error::Index {
            index: dists.iter().copied().max().unwrap_or_default(),
            len: self.bins.len() as u32,
        };
        let bin = self.find_bin_idx(&dists, scaffolds).ok_or_else(not_found)?;
        Ok(start_idx + offset + bin as u64)
    }

    /// The bit's number of points, its features (indices into
    /// [`get_feat_families`](SigFactory::get_feat_families)), and its
    /// scaffold (a bin index for each distance). `Error::Index` past the
    /// signature. Mirrors `SigFactory.GetBitInfo`.
    pub fn get_bit_info(&self, idx: u64) -> Result<(usize, Vec<usize>, Vec<usize>)> {
        if idx >= self.sig_size {
            return Err(Error::Index {
                index: u32::try_from(idx).unwrap_or(u32::MAX),
                len: u32::try_from(self.sig_size).unwrap_or(u32::MAX),
            });
        }
        let start = |n: usize| -> Result<u64> {
            self.starts.get(&n).copied().ok_or_else(|| Error::Value {
                message: "the factory needs init".to_owned(),
            })
        };
        let mut n_pts = self.min_point_count;
        while n_pts < self.max_point_count && start(n_pts + 1)? <= idx {
            n_pts += 1;
        }
        let offset_from_start = idx - start(n_pts)?;
        let n_dists = point_dists(n_pts)?.len();
        let scaffolds = self.scaffolds.get(&n_dists).ok_or_else(|| Error::Value {
            message: "the factory needs init".to_owned(),
        })?;
        let n_scaffolds = scaffolds.len() as u64;
        let proto_idx = (offset_from_start / n_scaffolds) as usize;
        let combo = get_index_combinations(self.n_feats, n_pts)
            .into_iter()
            .nth(proto_idx)
            .ok_or_else(|| Error::Value {
                message: "the factory needs init".to_owned(),
            })?;
        let scaffold = scaffolds[(offset_from_start % n_scaffolds) as usize].clone();
        Ok((n_pts, combo, scaffold))
    }

    /// Computes the signature's layout from the parameters; call it after
    /// changing them. `Error::Key` for a `max_point_count` of 10 or more
    /// (Python's dictionary has no distances past 10 points). Mirrors
    /// `SigFactory.Init`.
    pub fn init(&mut self) -> Result<()> {
        point_dists(self.max_point_count + 1)?;
        self.scaffolds.clear();
        self.starts.clear();
        self.n_feats = self
            .feat_factory
            .get_feature_families()
            .iter()
            .filter(|fam| !self.skip_feats.contains(fam))
            .count();
        let mut accum = 0u64;
        for i in self.min_point_count..=self.max_point_count {
            self.starts.insert(i, accum);
            let n_dists_here = point_dists(i)?.len();
            let scaffolds_here = get_possible_scaffolds_with_use_triangle_inequality(
                i,
                &self.bins,
                self.triangle_prune_bins,
            )?;
            accum += num_combinations(self.n_feats, i) * scaffolds_here.len() as u64;
            self.scaffolds.insert(n_dists_here, scaffolds_here);
        }
        self.sig_size = accum;
        Ok(())
    }

    /// The number of bits in a signature. Mirrors `SigFactory.GetSigSize`.
    pub fn get_sig_size(&self) -> u64 {
        self.sig_size
    }
}
