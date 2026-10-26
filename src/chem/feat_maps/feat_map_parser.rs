// Copyright (c) 2026 Steven Kearnes
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

//! Ports `rdkit.Chem.FeatMaps.FeatMapParser`, a Python module: the
//! text format of feature maps.
//!
//! ```text
//! ScoreMode=All
//! DirScoreMode=Ignore
//!
//! BeginParams
//!   family=Aromatic radius=2.5 width=1.0 profile=Gaussian
//!   family=Acceptor radius=1.5
//! EndParams
//!
//! # optional
//! BeginPoints
//!   family=Acceptor pos=(1.0, 0.0, 5.0) weight=1.25 dir=(1, 1, 0)
//!   family=Aromatic pos=(0.0,1.0,0.0) weight=2.0 dir=(0,0,1) dir=(0,0,-1)
//!   family=Acceptor pos=(1.0,1.0,2.0) weight=1.25
//! EndPoints
//! ```
//!
//! `#` starts a comment, and a line ending in `\` continues on the next.
//! Python's `FeatMapParseError` is `Error::Parse` here.

use std::collections::BTreeMap;
use std::io::Read;

use crate::chem::feat_maps::feat_map_point::FeatMapPoint;
use crate::chem::feat_maps::feat_maps::{
    FeatDirScoreMode, FeatMap, FeatMapParams, FeatMapScoreMode, FeatProfile,
};
use crate::geometry::rd_geometry::Point3D;
use crate::py::parse_float;
use crate::{Error, Result};

fn parse_error(message: String) -> Error {
    Error::Parse { message }
}

/// The end of a run of ASCII letters at `i`, or `i`.
fn letters_end(chars: &[char], i: usize) -> usize {
    let mut j = i;
    while j < chars.len() && chars[j].is_ascii_alphabetic() {
        j += 1;
    }
    j
}

/// The end of a run of spaces at `i`, or `i`.
fn spaces_end(chars: &[char], i: usize) -> usize {
    let mut j = i;
    while j < chars.len() && chars[j] == ' ' {
        j += 1;
    }
    j
}

/// A match of `([a-zA-Z]+) *= *` at `i`: the name's end and the match's.
fn name_eq_at(chars: &[char], i: usize) -> Option<(usize, usize)> {
    let name_end = letters_end(chars, i);
    if name_end == i {
        return None;
    }
    let eq = spaces_end(chars, name_end);
    if chars.get(eq) != Some(&'=') {
        return None;
    }
    Some((name_end, spaces_end(chars, eq + 1)))
}

/// Python's `re.findall(r'([a-zA-Z]+) *= *(\S+)', text)`.
fn find_name_values(text: &str) -> Vec<(String, String)> {
    let chars: Vec<char> = text.chars().collect();
    let mut res = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        if let Some((name_end, value_start)) = name_eq_at(&chars, i) {
            let mut value_end = value_start;
            while value_end < chars.len() && !chars[value_end].is_whitespace() {
                value_end += 1;
            }
            if value_end > value_start {
                res.push((
                    chars[i..name_end].iter().collect(),
                    chars[value_start..value_end].iter().collect(),
                ));
                i = value_end;
                continue;
            }
        }
        i += 1;
    }
    res
}

/// Python's `re.split(r'([a-zA-Z]+) *= *', text)` without its empty
/// strings: text between matches, and each match's name.
fn split_names(text: &str) -> Vec<String> {
    let chars: Vec<char> = text.chars().collect();
    let mut res = Vec::new();
    let (mut i, mut piece_start) = (0, 0);
    while i < chars.len() {
        if let Some((name_end, end)) = name_eq_at(&chars, i) {
            res.push(chars[piece_start..i].iter().collect());
            res.push(chars[i..name_end].iter().collect());
            piece_start = end;
            i = end;
        } else {
            i += 1;
        }
    }
    res.push(chars[piece_start..].iter().collect());
    res.retain(|s: &String| !s.is_empty());
    res
}

/// Reads feature maps from their text format. Mirrors
/// `FeatMapParser.FeatMapParser`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FeatMapParser {
    /// The lines to parse. Mirrors `FeatMapParser.data`.
    pub data: Vec<String>,
    line_num: usize,
}

impl FeatMapParser {
    /// A parser of a file's lines; `Error::Os` if it cannot be read and
    /// `Error::Value` if it is not UTF-8. Mirrors `FeatMapParser(file)`.
    pub fn from_file<R: Read + ?Sized>(file: &mut R) -> Result<FeatMapParser> {
        let mut bytes = Vec::new();
        file.read_to_end(&mut bytes).map_err(|err| Error::Os {
            message: err.to_string(),
        })?;
        let text = String::from_utf8(bytes).map_err(|err| Error::Value {
            message: format!("the feature map is not UTF-8: {err}"),
        })?;
        // Python's readlines keeps the line ends, which parsing strips.
        Ok(FeatMapParser {
            data: text.split_inclusive('\n').map(str::to_owned).collect(),
            line_num: 0,
        })
    }

    /// A parser of text. Mirrors `FeatMapParser(data=data)`.
    pub fn from_data(data: &str) -> FeatMapParser {
        let mut parser = FeatMapParser::default();
        parser.set_data(data);
        parser
    }

    /// Parses `data`'s lines from the start. Mirrors `FeatMapParser.SetData`.
    pub fn set_data(&mut self, data: &str) {
        self.set_data_lines(data.split('\n').map(str::to_owned).collect());
    }

    /// Parses these lines from the start. Mirrors `FeatMapParser.SetData`
    /// with a list of lines.
    pub fn set_data_lines(&mut self, lines: Vec<String>) {
        self.data = lines;
        self.line_num = 0;
    }

    /// The next line with content, without its comment, joined with the
    /// lines it continues onto.
    fn next_line(&mut self) -> String {
        let mut txt = String::new();
        while let Some(line) = self.data.get(self.line_num) {
            let l = line.split('#').next().unwrap_or_default().trim();
            self.line_num += 1;
            if !l.is_empty() {
                txt += l;
                if !l.ends_with('\\') {
                    break;
                }
            }
        }
        txt
    }

    /// Parses a feature map. Unknown keywords and modes are `Error::Parse`.
    /// Mirrors `FeatMapParser.Parse`.
    pub fn parse(&mut self) -> Result<FeatMap> {
        self.parse_with_feat_map(FeatMap::default())
    }

    /// Parses into `feat_map`, adding its points and replacing its modes and
    /// parameters, and returns it. Mirrors `FeatMapParser.Parse(featMap)`.
    pub fn parse_with_feat_map(&mut self, mut feat_map: FeatMap) -> Result<FeatMap> {
        let mut l = self.next_line().trim().to_owned();
        while !l.is_empty() {
            let split_l: Vec<&str> = l.split('=').collect();
            let keyword = split_l[0].trim().to_lowercase();
            if split_l.len() == 1 {
                match keyword.as_str() {
                    "beginpoints" => {
                        for pt in self.parse_feat_point_block()? {
                            feat_map.add_feat_point(pt)?;
                        }
                    }
                    "beginparams" => feat_map.params = self.parse_param_block()?,
                    _ => {
                        return Err(parse_error(format!(
                            "Unrecognized keyword {keyword} on line {}",
                            self.line_num
                        )));
                    }
                }
            } else {
                let val = split_l[1].trim();
                match keyword.as_str() {
                    "scoremode" => {
                        feat_map.score_mode = match val {
                            "All" => FeatMapScoreMode::All,
                            "Closest" => FeatMapScoreMode::Closest,
                            "Best" => FeatMapScoreMode::Best,
                            _ => {
                                return Err(parse_error(format!(
                                    "ScoreMode {val} not recognized on line {}",
                                    self.line_num
                                )));
                            }
                        }
                    }
                    "dirscoremode" => {
                        feat_map.dir_score_mode = match val {
                            "Ignore" => FeatDirScoreMode::Ignore,
                            "DotFullRange" => FeatDirScoreMode::DotFullRange,
                            "DotPosRange" => FeatDirScoreMode::DotPosRange,
                            _ => {
                                return Err(parse_error(format!(
                                    "DirScoreMode {val} not recognized on line {}",
                                    self.line_num
                                )));
                            }
                        }
                    }
                    _ => {
                        return Err(parse_error(format!(
                            "Unrecognized keyword {keyword} on line {}",
                            self.line_num
                        )));
                    }
                }
            }
            l = self.next_line().trim().to_owned();
        }
        Ok(feat_map)
    }

    /// Parses the lines up to `EndParams` into scoring parameters by
    /// family; a line without a family takes the last line's, and a block
    /// without its end is `Error::Parse`. Mirrors
    /// `FeatMapParser.ParseParamBlock`.
    pub fn parse_param_block(&mut self) -> Result<BTreeMap<String, FeatMapParams>> {
        let mut params = BTreeMap::new();
        let mut family: Option<String> = None;
        let mut l = self.next_line();
        while !l.is_empty() && l != "EndParams" {
            let mut param = FeatMapParams::default();
            for (name, val) in find_name_values(&l) {
                match name.to_lowercase().as_str() {
                    "family" => family = Some(val),
                    "radius" => param.radius = parse_float(&val)?,
                    "width" => param.width = parse_float(&val)?,
                    "profile" => {
                        param.feat_profile = match val.as_str() {
                            "Gaussian" => FeatProfile::Gaussian,
                            "Triangle" => FeatProfile::Triangle,
                            "Box" => FeatProfile::Box,
                            _ => {
                                return Err(parse_error(format!(
                                    "Profile {val} not recognized on line {}",
                                    self.line_num
                                )));
                            }
                        }
                    }
                    other => {
                        return Err(parse_error(format!(
                            "FeatMapParam option {other} not recognized on line {}",
                            self.line_num
                        )));
                    }
                }
            }
            // Python reads the family before any line gives one as an
            // unbound local.
            let family = family
                .clone()
                .ok_or_else(|| parse_error(format!("no family on line {}", self.line_num)))?;
            params.insert(family, param);
            l = self.next_line();
        }
        if l != "EndParams" {
            return Err(parse_error("EndParams line not found".to_owned()));
        }
        Ok(params)
    }

    /// A point as `(x, y, z)`, with or without the parentheses; another
    /// number of coordinates is `Error::Value`.
    fn parse_point(txt: &str) -> Result<Point3D> {
        let mut txt = txt.trim();
        if let Some(rest) = txt.strip_prefix('(') {
            txt = rest;
        }
        if let Some(rest) = txt.strip_suffix(')') {
            txt = rest;
        }
        let split_l: Vec<&str> = txt.split(',').collect();
        if split_l.len() != 3 {
            return Err(Error::Value {
                message: "Bad location string".to_owned(),
            });
        }
        Ok(Point3D::new(
            parse_float(split_l[0])?,
            parse_float(split_l[1])?,
            parse_float(split_l[2])?,
        ))
    }

    /// Parses the lines up to `EndPoints`, or the end, into points.
    /// Unknown options are `Error::Parse`. Mirrors
    /// `FeatMapParser.ParseFeatPointBlock`.
    pub fn parse_feat_point_block(&mut self) -> Result<Vec<FeatMapPoint>> {
        let mut feats = Vec::new();
        let mut l = self.next_line();
        while !l.is_empty() && l != "EndPoints" {
            let vals = split_names(&l);
            let mut p = FeatMapPoint::default();
            for pair in vals.chunks(2) {
                let name = pair[0].to_lowercase();
                // Python indexes past the end for an option without a value.
                let value = pair.get(1).ok_or_else(|| Error::Value {
                    message: format!("no value for {name} on line {}", self.line_num),
                })?;
                match name.as_str() {
                    "family" => p.set_family(value.trim()),
                    "weight" => p.weight = parse_float(value)?,
                    "pos" => p.set_pos(Self::parse_point(value)?),
                    "dir" => p.feat_dirs.push(Self::parse_point(value)?),
                    _ => {
                        return Err(parse_error(format!(
                            "FeatPoint option {name} not recognized on line {}",
                            self.line_num
                        )));
                    }
                }
            }
            feats.push(p);
            l = self.next_line();
        }
        Ok(feats)
    }
}
