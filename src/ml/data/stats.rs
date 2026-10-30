// Copyright (c) 2026 Genesis Molecular AI, Inc.
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.ML.Data.Stats`, a Python module: statistics on data
//! matrices (rows of points).
//!
//! Sums follow Python's (its built-in `sum`, row by row), so the means,
//! deviations, and correlations match Python's to the last bit. Python's
//! matrix products go through BLAS and its eigenvectors through LAPACK, whose
//! rounding and eigenvector signs vary by platform; here they are plain
//! loops and a Jacobi eigensolver, so those results match to within
//! rounding, and eigenvectors up to sign.

use crate::{Error, Result};

/// Python's `sum(rows)`: the rows added item by item, from zero.
fn column_sums(mat: &[Vec<f64>]) -> Vec<f64> {
    let ncols = mat.first().map_or(0, Vec::len);
    let mut sums = vec![0.0; ncols];
    for row in mat {
        for (s, &x) in sums.iter_mut().zip(row) {
            *s += x;
        }
    }
    sums
}

/// Python's `sum(values)`.
fn seq_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    values.into_iter().fold(0.0, |acc, x| acc + x)
}

fn check_rectangular(mat: &[Vec<f64>]) -> Result<usize> {
    let ncols = mat.first().map_or(0, Vec::len);
    if mat.iter().any(|row| row.len() != ncols) {
        return Err(Error::Value {
            message: "the matrix's rows have different lengths".to_owned(),
        });
    }
    Ok(ncols)
}

/// The matrix with each column's mean subtracted and the result divided by
/// the column's sample standard deviation. Mirrors
/// `Stats.StandardizeMatrix`, which also changes the matrix it is given.
///
/// In RDKit 2025.09.6 with NumPy 2, Python's raises `TypeError` for any
/// two-dimensional matrix (it takes `math.sqrt` of the deviations' array)
/// and standardizes only one-dimensional vectors; this standardizes each
/// column as Python does such a vector. A column with no deviation divides
/// by zero, to infinities and NaNs, as NumPy does.
///
/// `Error::Value` for rows of different lengths.
pub fn standardize_matrix(mat: &[Vec<f64>]) -> Result<Vec<Vec<f64>>> {
    check_rectangular(mat)?;
    let n_objs = mat.len() as f64;
    let avgs: Vec<f64> = column_sums(mat).into_iter().map(|s| s / n_objs).collect();
    let centered: Vec<Vec<f64>> = mat
        .iter()
        .map(|row| row.iter().zip(&avgs).map(|(x, a)| x - a).collect())
        .collect();
    let squares: Vec<Vec<f64>> = centered
        .iter()
        .map(|row| row.iter().map(|x| x * x).collect())
        .collect();
    let devs: Vec<f64> = column_sums(&squares)
        .into_iter()
        .map(|s| (s / (n_objs - 1.0)).sqrt())
        .collect();
    Ok(centered
        .into_iter()
        .map(|row| row.into_iter().zip(&devs).map(|(x, d)| x / d).collect())
        .collect())
}

/// The sample covariance matrix of the columns. Mirrors
/// `Stats.FormCovarianceMatrix`, which also centers the matrix it is given.
///
/// `Error::Value` for rows of different lengths.
pub fn form_covariance_matrix(mat: &[Vec<f64>]) -> Result<Vec<Vec<f64>>> {
    let ncols = check_rectangular(mat)?;
    let n_pts = mat.len() as f64;
    let means: Vec<f64> = column_sums(mat).into_iter().map(|s| s / n_pts).collect();
    let centered: Vec<Vec<f64>> = mat
        .iter()
        .map(|row| row.iter().zip(&means).map(|(x, m)| x - m).collect())
        .collect();
    Ok((0..ncols)
        .map(|i| {
            (0..ncols)
                .map(|j| seq_sum(centered.iter().map(|row| row[i] * row[j])) / (n_pts - 1.0))
                .collect()
        })
        .collect())
}

/// The matrix of Pearson correlations between the columns, with 0 where a
/// column does not vary. Mirrors `Stats.FormCorrelationMatrix`.
///
/// `Error::Value` for rows of different lengths.
#[allow(clippy::needless_range_loop)]
pub fn form_correlation_matrix(mat: &[Vec<f64>]) -> Result<Vec<Vec<f64>>> {
    let n_vars = check_rectangular(mat)?;
    let n = mat.len() as f64;
    let column = |i: usize| mat.iter().map(move |row| row[i]);
    let mut res = vec![vec![0.0; n_vars]; n_vars];
    for i in 0..n_vars {
        let sum_x = seq_sum(column(i));
        let sum_x2 = seq_sum(column(i).map(|x| x * x));
        for j in i..n_vars {
            let sum_y = seq_sum(column(j));
            let sum_y2 = seq_sum(column(j).map(|y| y * y));
            let sum_xy = seq_sum(column(i).zip(column(j)).map(|(x, y)| x * y));
            let numerator = n * sum_xy - sum_x * sum_y;
            let denom = ((n * sum_x2 - sum_x * sum_x) * (n * sum_y2 - sum_y * sum_y)).sqrt();
            let r = if denom != 0.0 { numerator / denom } else { 0.0 };
            res[i][j] = r;
            res[j][i] = r;
        }
    }
    Ok(res)
}

/// The eigenvalues and eigenvectors (as columns) of a symmetric matrix, by
/// cyclic Jacobi rotations. Each eigenvector has unit length and its
/// largest component (the first, among equals) positive.
#[allow(clippy::needless_range_loop)]
pub(crate) fn symmetric_eigen(mat: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = mat.len();
    let mut a: Vec<Vec<f64>> = mat.to_vec();
    let mut v: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    for _sweep in 0..100 {
        let off: f64 = (0..n)
            .flat_map(|i| (0..n).filter(move |&j| j != i).map(move |j| (i, j)))
            .map(|(i, j)| a[i][j] * a[i][j])
            .sum();
        if off == 0.0 {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                if a[p][q] == 0.0 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for row in v.iter_mut() {
                    let (vkp, vkq) = (row[p], row[q]);
                    row[p] = c * vkp - s * vkq;
                    row[q] = s * vkp + c * vkq;
                }
            }
        }
    }
    for j in 0..n {
        let mut largest = 0;
        for i in 0..n {
            if v[i][j].abs() > v[largest][j].abs() {
                largest = i;
            }
        }
        if v[largest][j] < 0.0 {
            for row in v.iter_mut() {
                row[j] = -row[j];
            }
        }
    }
    ((0..n).map(|i| a[i][i]).collect(), v)
}

/// A principal components analysis of the columns: the eigenvalues of their
/// correlation matrix, largest first, and the matching eigenvectors, as
/// rows. Mirrors `Stats.PrincipalComponents(mat)`.
///
/// Python's takes rows, not columns, of LAPACK's matrix of eigenvectors,
/// picking them by each eigenvalue's place in LAPACK's order, so its
/// "eigenvectors" depend on that order and are not eigenvectors; this
/// returns the eigenvectors. Each is signed so that its largest component
/// is positive, where LAPACK's signs are arbitrary.
///
/// `Error::Value` for rows of different lengths.
pub fn principal_components(mat: &[Vec<f64>]) -> Result<(Vec<f64>, Vec<Vec<f64>>)> {
    principal_components_with_reverse_order(mat, true)
}

/// As [`principal_components`], smallest eigenvalue first unless
/// `reverse_order`. Mirrors `Stats.PrincipalComponents(mat, reverseOrder)`.
pub fn principal_components_with_reverse_order(
    mat: &[Vec<f64>],
    reverse_order: bool,
) -> Result<(Vec<f64>, Vec<Vec<f64>>)> {
    let cov = form_correlation_matrix(mat)?;
    let (vals, vects) = symmetric_eigen(&cov);
    let mut order: Vec<usize> = (0..vals.len()).collect();
    order.sort_by(|&a, &b| vals[a].total_cmp(&vals[b]));
    if reverse_order {
        order.reverse();
    }
    Ok((
        order.iter().map(|&x| vals[x]).collect(),
        order
            .iter()
            .map(|&x| vects.iter().map(|row| row[x]).collect())
            .collect(),
    ))
}

/// The points, centered on their mean, multiplied by `t_form_mat`. Mirrors
/// `Stats.TransformPoints`.
///
/// `Error::Value` for points of different lengths, or a matrix whose rows
/// are not as long as the points.
pub fn transform_points(t_form_mat: &[Vec<f64>], pts: &[Vec<f64>]) -> Result<Vec<Vec<f64>>> {
    let dim = check_rectangular(pts)?;
    if t_form_mat.iter().any(|row| row.len() != dim) {
        return Err(Error::Value {
            message: "the matrix's rows must be as long as the points".to_owned(),
        });
    }
    let n_pts = pts.len() as f64;
    let avg: Vec<f64> = column_sums(pts).into_iter().map(|s| s / n_pts).collect();
    Ok(pts
        .iter()
        .map(|p| {
            let centered: Vec<f64> = p.iter().zip(&avg).map(|(x, a)| x - a).collect();
            t_form_mat
                .iter()
                .map(|row| seq_sum(row.iter().zip(&centered).map(|(m, x)| m * x)))
                .collect()
        })
        .collect())
}

/// The mean and sample standard deviation of the values, `(0, 0)` for
/// none, and a deviation of 0 for one. Mirrors `Stats.MeanAndDev(vect)`.
pub fn mean_and_dev(vect: &[f64]) -> (f64, f64) {
    mean_and_dev_with_sample_sd(vect, true)
}

/// As [`mean_and_dev`], with the population standard deviation unless
/// `sample_sd`. Mirrors `Stats.MeanAndDev(vect, sampleSD)`.
pub fn mean_and_dev_with_sample_sd(vect: &[f64], sample_sd: bool) -> (f64, f64) {
    let n = vect.len();
    if n == 0 {
        return (0.0, 0.0);
    }
    let mean = seq_sum(vect.iter().copied()) / n as f64;
    let ss = seq_sum(vect.iter().map(|x| (x - mean) * (x - mean)));
    let dev = if n > 1 {
        if sample_sd {
            (ss / (n - 1) as f64).sqrt()
        } else {
            (ss / n as f64).sqrt()
        }
    } else {
        0.0
    };
    (mean, dev)
}

/// The coefficient of determination of predictions with the residual sum
/// of squares `resid_sum`: 1 minus it over the values' sum of squared
/// deviations from their mean. Mirrors `Stats.R2`, which returns the tuple
/// `(0, 0)` for no values; this returns 0 there.
pub fn r2(orig: &[f64], resid_sum: f64) -> f64 {
    let n = orig.len();
    if n == 0 {
        return 0.0;
    }
    let o_mean = seq_sum(orig.iter().copied()) / n as f64;
    let o_var = seq_sum(orig.iter().map(|x| (x - o_mean) * (x - o_mean)));
    1.0 - resid_sum / o_var
}

// The table's 2.718 is not e.
#[allow(clippy::approx_constant)]
const T_TABLE: [[f64; 8]; 48] = [
    [1.0, 3.078, 6.314, 12.71, 31.82, 63.66, 318.30, 637.0],
    [2.0, 1.886, 2.920, 4.303, 6.965, 9.925, 22.330, 31.6],
    [3.0, 1.638, 2.353, 3.182, 4.541, 5.841, 10.210, 12.92],
    [4.0, 1.533, 2.132, 2.776, 3.747, 4.604, 7.173, 8.610],
    [5.0, 1.476, 2.015, 2.571, 3.365, 4.032, 5.893, 6.869],
    [6.0, 1.440, 1.943, 2.447, 3.143, 3.707, 5.208, 5.959],
    [7.0, 1.415, 1.895, 2.365, 2.998, 3.499, 4.785, 5.408],
    [8.0, 1.397, 1.860, 2.306, 2.896, 3.355, 4.501, 5.041],
    [9.0, 1.383, 1.833, 2.262, 2.821, 3.250, 4.297, 4.781],
    [10.0, 1.372, 1.812, 2.228, 2.764, 3.169, 4.144, 4.587],
    [11.0, 1.363, 1.796, 2.201, 2.718, 3.106, 4.025, 4.437],
    [12.0, 1.356, 1.782, 2.179, 2.681, 3.055, 3.930, 4.318],
    [13.0, 1.350, 1.771, 2.160, 2.650, 3.012, 3.852, 4.221],
    [14.0, 1.345, 1.761, 2.145, 2.624, 2.977, 3.787, 4.140],
    [15.0, 1.341, 1.753, 2.131, 2.602, 2.947, 3.733, 4.073],
    [16.0, 1.337, 1.746, 2.120, 2.583, 2.921, 3.686, 4.015],
    [17.0, 1.333, 1.740, 2.110, 2.567, 2.898, 3.646, 3.965],
    [18.0, 1.330, 1.734, 2.101, 2.552, 2.878, 3.610, 3.922],
    [19.0, 1.328, 1.729, 2.093, 2.539, 2.861, 3.579, 3.883],
    [20.0, 1.325, 1.725, 2.086, 2.528, 2.845, 3.552, 3.850],
    [21.0, 1.323, 1.721, 2.080, 2.518, 2.831, 3.527, 3.819],
    [22.0, 1.321, 1.717, 2.074, 2.508, 2.819, 3.505, 3.792],
    [23.0, 1.319, 1.714, 2.069, 2.500, 2.807, 3.485, 3.768],
    [24.0, 1.318, 1.711, 2.064, 2.492, 2.797, 3.467, 3.745],
    [25.0, 1.316, 1.708, 2.060, 2.485, 2.787, 3.450, 3.725],
    [26.0, 1.315, 1.706, 2.056, 2.479, 2.779, 3.435, 3.707],
    [27.0, 1.314, 1.703, 2.052, 2.473, 2.771, 3.421, 3.690],
    [28.0, 1.313, 1.701, 2.048, 2.467, 2.763, 3.408, 3.674],
    [29.0, 1.311, 1.699, 2.045, 2.462, 2.756, 3.396, 3.659],
    [30.0, 1.310, 1.697, 2.042, 2.457, 2.750, 3.385, 3.646],
    [32.0, 1.309, 1.694, 2.037, 2.449, 2.738, 3.365, 3.622],
    [34.0, 1.307, 1.691, 2.032, 2.441, 2.728, 3.348, 3.601],
    [36.0, 1.306, 1.688, 2.028, 2.434, 2.719, 3.333, 3.582],
    [38.0, 1.304, 1.686, 2.024, 2.429, 2.712, 3.319, 3.566],
    [40.0, 1.303, 1.684, 2.021, 2.423, 2.704, 3.307, 3.551],
    [42.0, 1.302, 1.682, 2.018, 2.418, 2.698, 3.296, 3.538],
    [44.0, 1.301, 1.680, 2.015, 2.414, 2.692, 3.286, 3.526],
    [46.0, 1.300, 1.679, 2.013, 2.410, 2.687, 3.277, 3.515],
    [48.0, 1.299, 1.677, 2.011, 2.407, 2.682, 3.269, 3.505],
    [50.0, 1.299, 1.676, 2.009, 2.403, 2.678, 3.261, 3.496],
    [55.0, 1.297, 1.673, 2.004, 2.396, 2.668, 3.245, 3.476],
    [60.0, 1.296, 1.671, 2.000, 2.390, 2.660, 3.232, 3.460],
    [65.0, 1.295, 1.669, 1.997, 2.385, 2.654, 3.220, 3.447],
    [70.0, 1.294, 1.667, 1.994, 2.381, 2.648, 3.211, 3.435],
    [80.0, 1.292, 1.664, 1.990, 2.374, 2.639, 3.195, 3.416],
    [100.0, 1.290, 1.660, 1.984, 2.364, 2.626, 3.174, 3.390],
    [150.0, 1.287, 1.655, 1.976, 2.351, 2.609, 3.145, 3.357],
    [200.0, 1.286, 1.653, 1.972, 2.345, 2.601, 3.131, 3.340],
];

/// The two-sided Student's t table: each row is the degrees of freedom,
/// then t at each confidence level of [`t_confs`]. Mirrors `Stats.tTable`.
pub fn t_table() -> &'static [[f64; 8]] {
    &T_TABLE
}

/// The confidence levels (percent) [`t_table`] covers, each with its
/// column. Mirrors `Stats.tConfs`.
pub fn t_confs() -> Vec<(f64, usize)> {
    vec![
        (80.0, 1),
        (90.0, 2),
        (95.0, 3),
        (98.0, 4),
        (99.0, 5),
        (99.8, 6),
        (99.9, 7),
    ]
}

/// The half-width of the 95% confidence interval of a mean of `n` values
/// with standard deviation `sd`, from [`t_table`]'s first row with at least
/// `n - 1` degrees of freedom (its last row past the table). Mirrors
/// `Stats.GetConfidenceInterval(sd, n)`.
pub fn get_confidence_interval(sd: f64, n: i64) -> Result<f64> {
    get_confidence_interval_with_level(sd, n, 95.0)
}

/// As [`get_confidence_interval`], at the confidence `level`, one of
/// [`t_confs`]' levels. Mirrors `Stats.GetConfidenceInterval(sd, n,
/// level)`.
///
/// `Error::Key` for another level, as Python raises.
pub fn get_confidence_interval_with_level(sd: f64, n: i64, level: f64) -> Result<f64> {
    let col = t_confs()
        .into_iter()
        .find(|&(l, _)| l == level)
        .map(|(_, col)| col)
        .ok_or_else(|| Error::Key {
            key: level.to_string(),
        })?;
    let dofs = (n - 1) as f64;
    let sem = sd / (n as f64).sqrt();
    let table = t_table();
    let row = table
        .iter()
        .find(|row| row[0].partial_cmp(&dofs) != Some(std::cmp::Ordering::Less))
        .unwrap_or(&table[table.len() - 1]);
    Ok(row[col] * sem)
}
