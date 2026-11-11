//! PWM scanning with an exact score distribution (discretized) for p-values.

use std::collections::HashMap;

use crate::promoter::reverse_complement;
use crate::pwm::Pwm;

#[derive(Clone, Debug)]
pub struct RawHit {
    /// 1-based inclusive position in the oriented promoter sequence.
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub score: f64,
    pub p_value: f64,
}

#[derive(Clone, Debug)]
pub struct PwmNull {
    /// Sorted (discretized score, cumulative probability up to this score).
    sorted: Vec<(i32, f64)>,
    scale: f64,
}

fn base_index(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

pub fn build_null(pwm: &Pwm, background: &[f64; 4], scale: f64) -> PwmNull {
    let log_odds = pwm.log_odds(background);
    let mut scores: Vec<[i32; 4]> = Vec::with_capacity(log_odds.len());
    for row in &log_odds {
        let mut r = [0_i32; 4];
        for b in 0..4 {
            r[b] = (row[b] * scale).round() as i32;
        }
        scores.push(r);
    }

    let mut dist: HashMap<i32, f64> = HashMap::new();
    dist.insert(0, 1.0);
    for row in &scores {
        let mut next: HashMap<i32, f64> = HashMap::with_capacity(dist.len() * 4);
        for (&s, &p) in &dist {
            for b in 0..4 {
                let bg = background[b];
                if bg <= 0.0 {
                    continue;
                }
                *next.entry(s + row[b]).or_insert(0.0) += p * bg;
            }
        }
        dist = next;
    }

    let mut sorted: Vec<(i32, f64)> = dist.into_iter().collect();
    sorted.sort_by_key(|(s, _)| *s);
    let mut cum = 0.0;
    for (_, p) in &mut sorted {
        cum += *p;
        *p = cum;
    }
    PwmNull { sorted, scale }
}

impl PwmNull {
    pub fn pvalue(&self, observed: f64) -> f64 {
        if self.sorted.is_empty() {
            return 1.0;
        }
        let obs = (observed * self.scale).round() as i32;
        let idx = self.sorted.partition_point(|(s, _)| *s < obs);
        if idx == 0 {
            return 1.0;
        }
        if idx >= self.sorted.len() {
            return 0.0;
        }
        let cum_before = self.sorted[idx - 1].1;
        (1.0 - cum_before).clamp(0.0, 1.0)
    }
}

fn scan_one_strand(
    seq: &[u8],
    log_odds: &[[f64; 4]],
    null: &PwmNull,
    p_threshold: f64,
    min_score: Option<f64>,
    strand: char,
) -> Vec<RawHit> {
    let l = log_odds.len();
    if seq.len() < l {
        return Vec::new();
    }
    let mut out = Vec::new();
    for offset in 0..=seq.len() - l {
        let mut score = 0.0;
        let mut valid = true;
        for (k, row) in log_odds.iter().enumerate() {
            let Some(b) = base_index(seq[offset + k]) else {
                valid = false;
                break;
            };
            score += row[b];
        }
        if !valid {
            continue;
        }
        if min_score.is_some_and(|m| score < m) {
            continue;
        }
        let p = null.pvalue(score);
        if p <= p_threshold {
            out.push(RawHit {
                start: offset as u64 + 1,
                end: (offset + l) as u64,
                strand,
                score,
                p_value: p,
            });
        }
    }
    out
}

/// Scan an oriented promoter on both strands. Hit coordinates refer to the
/// oriented promoter sequence (5' -> 3' in the gene's direction).
pub fn scan_promoter(
    promoter: &[u8],
    log_odds: &[[f64; 4]],
    null: &PwmNull,
    p_threshold: f64,
    min_score: Option<f64>,
    both_strands: bool,
) -> Vec<RawHit> {
    let mut hits = scan_one_strand(promoter, log_odds, null, p_threshold, min_score, '+');
    if both_strands {
        let n = promoter.len();
        let l = log_odds.len();
        let rc = reverse_complement(promoter);
        for hit in scan_one_strand(&rc, log_odds, null, p_threshold, min_score, '-') {
            let rc_start0 = (hit.start - 1) as usize;
            let start0 = n.saturating_sub(rc_start0 + l);
            hits.push(RawHit {
                start: start0 as u64 + 1,
                end: (start0 + l) as u64,
                strand: '-',
                score: hit.score,
                p_value: hit.p_value,
            });
        }
    }
    hits.sort_by(|a, b| {
        a.start
            .cmp(&b.start)
            .then(a.end.cmp(&b.end))
            .then(a.strand.cmp(&b.strand))
    });
    hits
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pwm() -> Pwm {
        Pwm {
            id: "M1".into(),
            name: "test".into(),
            matrix: vec![
                [0.97, 0.01, 0.01, 0.01],
                [0.01, 0.97, 0.01, 0.01],
            ],
        }
    }

    #[test]
    fn detects_consensus_and_pvalue_is_small() {
        let bg = [0.25; 4];
        let pwm = pwm();
        let log_odds = pwm.log_odds(&bg);
        let null = build_null(&pwm, &bg, 10.0);
        let hits = scan_promoter(b"TTACGTT", &log_odds, &null, 0.1, None, true);
        assert!(hits.iter().any(|h| h.strand == '+' && h.start == 3));
        for h in &hits {
            assert!(h.p_value <= 0.1);
        }
    }

    #[test]
    fn pvalue_is_monotonic() {
        let bg = [0.25; 4];
        let pwm = pwm();
        let null = build_null(&pwm, &bg, 10.0);
        let low = null.pvalue(0.0);
        let high = null.pvalue(20.0);
        assert!(low >= high);
    }
}
