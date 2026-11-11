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

/// Normalize a count column with a 0.5 pseudocount while avoiding overflow.
/// For columns whose maximum is small enough, direct arithmetic is safe; for
/// very large counts, scale by the column maximum first.
fn normalize_counts(col: [f64; 4]) -> [f64; 4] {
    let max = col.iter().copied().fold(0.0_f64, f64::max);
    if max <= 2.0 {
        let sum: f64 = col.iter().sum();
        let denom = sum + 2.0;
        col.map(|v| (v + 0.5) / denom)
    } else {
        let scaled = col.map(|v| v / max);
        let scaled_sum: f64 = scaled.iter().sum();
        let denom = scaled_sum + 2.0 / max;
        scaled.map(|v| (v + 0.5 / max) / denom)
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
    if rows
        .iter()
        .flatten()
        .any(|v| !v.is_finite() || *v < 0.0)
    {
        bail!("motif {id}: values must be finite and non-negative");
    }
    let max = rows.iter().flatten().copied().fold(0.0_f64, f64::max);
    let column_sums: Vec<f64> = (0..len)
        .map(|j| (0..4).map(|b| rows[b][j]).sum::<f64>())
        .collect();
    let looks_like_probabilities = max <= 1.0 + 1e-9
        && column_sums
            .iter()
            .all(|s| (*s - 1.0).abs() < 0.15);
    let mut matrix = Vec::with_capacity(len);
    for j in 0..len {
        let mut col = [0.0_f64; 4];
        if looks_like_probabilities {
            // Normalize by the actual column sum: JASPAR probability matrices
            // can be rounded and not sum exactly to 1.
            let denom = column_sums[j] + 4e-3;
            for (b, value) in col.iter_mut().enumerate() {
                *value = (rows[b][j] + 1e-3) / denom;
            }
        } else {
            col = normalize_counts([rows[0][j], rows[1][j], rows[2][j], rows[3][j]]);
        }
        matrix.push(col);
    }
    for (j, col) in matrix.iter().enumerate() {
        let sum: f64 = col.iter().sum();
        if col.iter().any(|v| !v.is_finite() || *v < 0.0)
            || !sum.is_finite()
            || (sum - 1.0).abs() > 1e-9
        {
            bail!("motif {id}: column {} failed to normalize (sum {sum})", j + 1);
        }
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
            let motif = id.as_deref().unwrap_or("<unnamed>");
            let v = token.parse::<f64>().with_context(|| {
                format!(
                    "motif {motif} row {} has invalid value '{token}' at line {}",
                    first as char,
                    lineno + 1
                )
            })?;
            values.push(v);
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

    #[test]
    fn normalizes_probability_columns_that_do_not_sum_to_one() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_pwm_norm_{}.txt",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, ">M1 TEST").unwrap();
            writeln!(f, "A [ 0.8 0.8 ]").unwrap();
            writeln!(f, "C [ 0.1 0.1 ]").unwrap();
            writeln!(f, "G [ 0.0 0.0 ]").unwrap();
            writeln!(f, "T [ 0.0 0.0 ]").unwrap();
        }
        let pwms = read_jaspar(&path).unwrap();
        for j in 0..2 {
            let s: f64 = pwms[0].matrix[j].iter().sum();
            assert!((s - 1.0).abs() < 1e-9, "column {j} sums to {s}");
        }
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn rejects_malformed_motif_tokens() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_pwm_malformed_{}.txt",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, ">M TEST").unwrap();
            writeln!(f, "A [ 10 oops 0 ]").unwrap();
            writeln!(f, "C [ 0 oops 10 ]").unwrap();
            writeln!(f, "G [ 0 oops 0 ]").unwrap();
            writeln!(f, "T [ 0 oops 0 ]").unwrap();
        }
        let err = format!("{:#}", read_jaspar(&path).unwrap_err());
        assert!(err.contains("invalid value 'oops'"), "{err}");
        assert!(err.contains("line 2"), "{err}");
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn huge_finite_counts_normalize_without_overflow() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_pwm_huge_{}.txt",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, ">M1 TEST").unwrap();
            writeln!(f, "A [ 1e308 ]").unwrap();
            writeln!(f, "C [ 1e308 ]").unwrap();
            writeln!(f, "G [ 1e308 ]").unwrap();
            writeln!(f, "T [ 1e308 ]").unwrap();
        }
        let pwms = read_jaspar(&path).unwrap();
        let col = pwms[0].matrix[0];
        let sum: f64 = col.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12, "sum={sum}");
        for v in col {
            assert!((v - 0.25).abs() < 1e-12, "v={v}");
        }
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn tiny_finite_counts_normalize_without_underflow() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_pwm_tiny_{}.txt",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, ">M1 TEST").unwrap();
            writeln!(f, "A [ 1e-308 ]").unwrap();
            writeln!(f, "C [ 1e-308 ]").unwrap();
            writeln!(f, "G [ 1e-308 ]").unwrap();
            writeln!(f, "T [ 1e-308 ]").unwrap();
        }
        let pwms = read_jaspar(&path).unwrap();
        let col = pwms[0].matrix[0];
        let sum: f64 = col.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12, "sum={sum}");
        for v in col {
            assert!((v - 0.25).abs() < 1e-12, "v={v}");
        }
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn rejects_negative_motif_values() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_pwm_negative_{}.txt",
            std::process::id()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, ">M1 TEST").unwrap();
            writeln!(f, "A [ 1.0 -0.1 ]").unwrap();
            writeln!(f, "C [ 0.0 0.1 ]").unwrap();
            writeln!(f, "G [ 0.0 0.0 ]").unwrap();
            writeln!(f, "T [ 0.0 0.0 ]").unwrap();
        }
        let err = read_jaspar(&path).unwrap_err().to_string();
        assert!(err.contains("non-negative"), "{err}");
        let _ = std::fs::remove_file(path);
    }
}
