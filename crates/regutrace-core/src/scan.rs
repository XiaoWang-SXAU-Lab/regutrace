//! PWM scanning with an exact score distribution (discretized) for p-values.

use std::collections::BTreeMap;

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
struct PwmNull {
    /// Sorted (discretized score, upper-tail probability P(null >= score)).
    /// The tail is accumulated directly in descending score order to avoid
    /// catastrophic cancellation for very small p-values.
    sorted: Vec<(i32, f64)>,
    /// Quantized per-position scores used to build `sorted`. Observed scores
    /// must be summed from this same table so the p-value tail is consistent.
    rows: Vec<[i32; 4]>,
}

/// A PWM bound to its log-odds table and null distribution. Constructing the
/// scanner from one PWM guarantees the two tables cannot be mismatched.
#[derive(Clone, Debug)]
pub struct Scanner {
    log_odds: Vec<[f64; 4]>,
    null: PwmNull,
}

impl Scanner {
    pub fn new(pwm: &Pwm, background: &[f64; 4], scale: f64) -> Self {
        Self {
            log_odds: pwm.log_odds(background),
            null: build_null(pwm, background, scale),
        }
    }

    pub fn scan(
        &self,
        promoter: &[u8],
        p_threshold: f64,
        min_score: Option<f64>,
        both_strands: bool,
    ) -> Vec<RawHit> {
        scan_promoter(
            promoter,
            &self.log_odds,
            &self.null,
            p_threshold,
            min_score,
            both_strands,
        )
    }
}

fn base_index(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn build_null(pwm: &Pwm, background: &[f64; 4], scale: f64) -> PwmNull {
    let log_odds = pwm.log_odds(background);
    let mut rows: Vec<[i32; 4]> = Vec::with_capacity(log_odds.len());
    for row in &log_odds {
        let mut r = [0_i32; 4];
        for b in 0..4 {
            r[b] = (row[b] * scale).round() as i32;
        }
        rows.push(r);
    }

    // BTreeMap makes the DP summation order deterministic across runs.
    let mut dist: BTreeMap<i32, f64> = BTreeMap::new();
    dist.insert(0, 1.0);
    for row in &rows {
        let mut next: BTreeMap<i32, f64> = BTreeMap::new();
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
    // Accumulate P(null >= score) directly from the largest score down. This
    // keeps tiny tail probabilities accurate instead of computing 1 - CDF.
    let mut tail = 0.0;
    for i in (0..sorted.len()).rev() {
        tail += sorted[i].1;
        sorted[i].1 = tail;
    }
    PwmNull { sorted, rows }
}

impl PwmNull {
    /// P(null discretized score >= `observed`). `observed` must be the sum of
    /// the same per-position quantized scores used to build the null.
    fn pvalue(&self, observed: i32) -> f64 {
        if self.sorted.is_empty() {
            return 1.0;
        }
        let idx = self.sorted.partition_point(|(s, _)| *s < observed);
        if idx >= self.sorted.len() {
            return 0.0;
        }
        self.sorted[idx].1.clamp(0.0, 1.0)
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
        let mut qscore = 0_i32;
        let mut valid = true;
        for (k, row) in log_odds.iter().enumerate() {
            let Some(b) = base_index(seq[offset + k]) else {
                valid = false;
                break;
            };
            score += row[b];
            qscore += null.rows[k][b];
        }
        if !valid {
            continue;
        }
        if min_score.is_some_and(|m| score < m) {
            continue;
        }
        let p = null.pvalue(qscore);
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
fn scan_promoter(
    promoter: &[u8],
    log_odds: &[[f64; 4]],
    null: &PwmNull,
    p_threshold: f64,
    min_score: Option<f64>,
    both_strands: bool,
) -> Vec<RawHit> {
    if log_odds.is_empty() || null.rows.len() != log_odds.len() {
        return Vec::new();
    }
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
    use std::collections::HashMap;

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
        let low = null.pvalue(0);
        let high = null.pvalue(200);
        assert!(low >= high);
    }

    #[test]
    fn null_distribution_matches_exact_enumeration() {
        let bg = [0.25; 4];
        let pwm = Pwm {
            id: "M2".into(),
            name: "test".into(),
            matrix: vec![
                [0.7, 0.1, 0.1, 0.1],
                [0.4, 0.3, 0.2, 0.1],
            ],
        };
        let scale = 10.0;
        let log_odds = pwm.log_odds(&bg);
        let null = build_null(&pwm, &bg, scale);
        let rows: Vec<[i32; 4]> = log_odds
            .iter()
            .map(|r| {
                let mut out = [0_i32; 4];
                for b in 0..4 {
                    out[b] = (r[b] * scale).round() as i32;
                }
                out
            })
            .collect();
        let mut dist: HashMap<i32, f64> = HashMap::new();
        for a in 0..4 {
            for b in 0..4 {
                let s = rows[0][a] + rows[1][b];
                *dist.entry(s).or_insert(0.0) += 0.25 * 0.25;
            }
        }
        for (&s, _) in &dist {
            let exact: f64 = dist
                .iter()
                .filter(|(t, _)| **t >= s)
                .map(|(_, &q)| q)
                .sum();
            let got = null.pvalue(s);
            assert!((got - exact).abs() < 1e-12, "score {s}: got {got}, exact {exact}");
        }
    }

    #[test]
    fn long_consensus_tail_keeps_precision() {
        let bg = [0.25; 4];
        let pwm = Pwm {
            id: "M5".into(),
            name: "A30".into(),
            matrix: vec![[0.97, 0.01, 0.01, 0.01]; 30],
        };
        let null = build_null(&pwm, &bg, 10.0);
        let q: i32 = (0..30).map(|k| null.rows[k][0]).sum();
        let p = null.pvalue(q);
        let exact = 0.25_f64.powi(30);
        assert!(
            (p - exact).abs() / exact < 1e-12,
            "p={p}, exact={exact}"
        );
    }

    #[test]
    fn scanner_binds_log_odds_and_null() {
        let bg = [0.25; 4];
        let scanner = Scanner::new(&pwm(), &bg, 10.0);
        let hits = scanner.scan(b"TTACGTT", 0.1, None, true);
        assert!(hits.iter().any(|h| h.strand == '+' && h.start == 3));
    }

    #[test]
    fn single_mismatch_pvalue_is_not_tail_collapsed() {
        let bg = [0.25; 4];
        let cons = 0.875;
        let misc = (1.0 - cons) / 3.0;
        // Motif ACGTAC with the consensus base at 87.5% and the other three at
        // 4.17%. A sequence with one mismatch should have p ~ 4.6e-3, not the
        // p of the perfect site (0.25^6 = 2.44e-4).
        let pwm = Pwm {
            id: "M4".into(),
            name: "ACGTAC".into(),
            matrix: vec![
                [cons, misc, misc, misc],
                [misc, cons, misc, misc],
                [misc, misc, cons, misc],
                [misc, misc, misc, cons],
                [cons, misc, misc, misc],
                [misc, cons, misc, misc],
            ],
        };
        let log_odds = pwm.log_odds(&bg);
        let null = build_null(&pwm, &bg, 10.0);
        // 5 consensus positions (one base at 25% background each) and one
        // mismatch (three bases, 75% total): P(>=64) = 6*0.25^5*0.75 + 0.25^6.
        let mut q = 0;
        for (k, &b) in b"ACGTAA".iter().enumerate() {
            let idx = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => unreachable!(),
            };
            q += null.rows[k][idx];
        }
        let p = null.pvalue(q);
        assert!((p - 0.004638671875).abs() < 1e-12, "p={p}");
        let hits = scan_promoter(b"ACGTAA", &log_odds, &null, 1e-3, None, true);
        assert!(hits.is_empty(), "{hits:?}");
    }

    #[test]
    fn scans_lowercase_sequence() {
        let bg = [0.25; 4];
        let pwm = pwm();
        let log_odds = pwm.log_odds(&bg);
        let null = build_null(&pwm, &bg, 10.0);
        let hits = scan_promoter(b"ttacgtt", &log_odds, &null, 0.1, None, true);
        assert!(hits.iter().any(|h| h.strand == '+' && h.start == 3));
    }

    #[test]
    fn minus_strand_hit_uses_oriented_promoter_coordinates() {
        let bg = [0.25; 4];
        let pwm = Pwm {
            id: "M3".into(),
            name: "AAGC".into(),
            matrix: vec![
                [0.97, 0.01, 0.01, 0.01],
                [0.97, 0.01, 0.01, 0.01],
                [0.01, 0.01, 0.97, 0.01],
                [0.01, 0.97, 0.01, 0.01],
            ],
        };
        let log_odds = pwm.log_odds(&bg);
        let null = build_null(&pwm, &bg, 10.0);
        // Forward promoter contains GCTT at 1-based 3..6 (0-based 2..5),
        // whose reverse complement is the motif AAGC.
        let hits = scan_promoter(b"TTGCTTTT", &log_odds, &null, 0.5, None, true);
        let minus = hits.iter().find(|h| h.strand == '-').expect("minus hit");
        assert_eq!(minus.start, 3);
        assert_eq!(minus.end, 6);
    }
}
