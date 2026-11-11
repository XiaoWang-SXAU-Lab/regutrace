//! JASPAR PFM/PWM reader and log-odds conversion.

use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::fasta::open_reader;

#[derive(Clone, Debug)]
pub struct Pwm {
    pub id: String,
    pub name: String,
    /// Positions x [A, C, G, T] probabilities.
    pub matrix: Vec<[f64; 4]>,
}

fn base_index(letter: u8) -> Option<usize> {
    match letter {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn finalize(id: String, name: String, rows: [Vec<f64>; 4]) -> Result<Pwm> {
    let len = rows[0].len();
    if len == 0 {
        bail!("motif {id} has no positions");
    }
    for (i, row) in rows.iter().enumerate() {
        if row.len() != len {
            bail!("motif {id}: row {} has {} positions, expected {len}", i, row.len());
        }
    }
    let max = rows.iter().flatten().copied().fold(0.0_f64, f64::max);
    let column_sums: Vec<f64> = (0..len)
        .map(|j| (0..4).map(|b| rows[b][j]).sum::<f64>())
        .collect();
    let looks_like_probabilities = max <= 1.0
        && column_sums
            .iter()
            .all(|s| (*s - 1.0).abs() < 0.15);
    let mut matrix = Vec::with_capacity(len);
    for j in 0..len {
        let mut col = [0.0_f64; 4];
        if looks_like_probabilities {
            for (b, value) in col.iter_mut().enumerate() {
                *value = (rows[b][j] + 1e-3) / (1.0 + 4e-3);
            }
        } else {
            let sum: f64 = (0..4).map(|b| rows[b][j]).sum();
            for (b, value) in col.iter_mut().enumerate() {
                *value = (rows[b][j] + 0.5) / (sum + 2.0);
            }
        }
        matrix.push(col);
    }
    Ok(Pwm { id, name, matrix })
}

/// Read a JASPAR-style PFM/PWM text file (multiple records allowed).
pub fn read_jaspar(path: &Path) -> Result<Vec<Pwm>> {
    let reader = BufReader::new(open_reader(path)?);
    let mut out = Vec::new();
    let mut id: Option<String> = None;
    let mut name = String::new();
    let mut rows: [Vec<f64>; 4] = Default::default();

    for (lineno, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("read motif line {}", lineno + 1))?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(prev_id) = id.take() {
                out.push(finalize(prev_id, std::mem::take(&mut name), std::mem::take(&mut rows))?);
            }
            let mut parts = rest.split_whitespace();
            let next_id = parts.next().unwrap_or("").to_string();
            if next_id.is_empty() {
                bail!("motif header without ID at line {}", lineno + 1);
            }
            id = Some(next_id);
            name = parts.collect::<Vec<_>>().join(" ");
            rows = Default::default();
            continue;
        }
        let Some(first) = line.as_bytes().first().copied() else {
            continue;
        };
        let Some(b) = base_index(first.to_ascii_uppercase()) else {
            continue;
        };
        let rest = &line[1..];
        let mut values = Vec::new();
        for token in rest.split(|c: char| c == '[' || c == ']' || c.is_whitespace()) {
            let token = token.trim();
            if token.is_empty() {
                continue;
            }
            if let Ok(v) = token.parse::<f64>() {
                values.push(v);
            }
        }
        if values.is_empty() {
            bail!("motif row {} has no numeric values at line {}", first as char, lineno + 1);
        }
        rows[b] = values;
    }
    if let Some(prev_id) = id.take() {
        out.push(finalize(prev_id, name, rows)?);
    }
    if out.is_empty() {
        bail!("no motifs found in {}", path.display());
    }
    Ok(out)
}

impl Pwm {
    pub fn log_odds(&self, background: &[f64; 4]) -> Vec<[f64; 4]> {
        self.matrix
            .iter()
            .map(|row| {
                let mut out = [0.0_f64; 4];
                for b in 0..4 {
                    let bg = background[b].max(1e-9);
                    out[b] = (row[b].max(1e-9) / bg).log2();
                }
                out
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn parses_counts_and_normalizes() {
        let path = std::env::temp_dir().join(format!("regutrace_pwm_{}.txt", std::process::id()));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, ">MA0001.1 TEST").unwrap();
            writeln!(f, "A [ 10 0 ]").unwrap();
            writeln!(f, "C [ 0 10 ]").unwrap();
            writeln!(f, "G [ 0 0 ]").unwrap();
            writeln!(f, "T [ 0 0 ]").unwrap();
        }
        let pwms = read_jaspar(&path).unwrap();
        assert_eq!(pwms.len(), 1);
        assert_eq!(pwms[0].id, "MA0001.1");
        assert_eq!(pwms[0].matrix.len(), 2);
        assert!(pwms[0].matrix[0][0] > 0.8);
        assert!(pwms[0].matrix[1][1] > 0.8);
        let _ = std::fs::remove_file(path);
    }
}
