//! Promoter extraction from a gene model and chromosome sequence.

use crate::gff::GeneModel;

#[derive(Clone, Debug)]
pub struct Promoter {
    pub gene_id: String,
    pub symbol: String,
    pub chrom: String,
    pub strand: char,
    /// 1-based inclusive genomic coordinates of the promoter window.
    pub start: u64,
    pub end: u64,
    /// Promoter sequence oriented 5' -> 3' in the gene's direction.
    pub seq: Vec<u8>,
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            _ => b'N',
        })
        .collect()
}

/// Extract `upstream_len` bases upstream of the representative TSS.
///
/// `chrom_seq` must be the full chromosome sequence (0-based, ACGTN).
pub fn extract_promoter(chrom_seq: &[u8], gene: &GeneModel, upstream_len: usize) -> Option<Promoter> {
    let (start0, end0, oriented) = if gene.strand == '-' {
        // tss is 1-based; for '-' it is the last base of the transcript, so the
        // 0-based index after the TSS base is `tss`.
        let start0 = gene.tss as usize;
        if start0 >= chrom_seq.len() {
            return None;
        }
        let end0 = start0.saturating_add(upstream_len).min(chrom_seq.len());
        if end0 <= start0 {
            return None;
        }
        (start0, end0, reverse_complement(&chrom_seq[start0..end0]))
    } else {
        // tss is 1-based; convert to the 0-based TSS base.
        let tss0 = gene.tss.saturating_sub(1) as usize;
        if tss0 == 0 || tss0 > chrom_seq.len() {
            return None;
        }
        let start0 = tss0.saturating_sub(upstream_len);
        let end0 = tss0;
        (start0, end0, chrom_seq[start0..end0].to_vec())
    };
    if oriented.is_empty() {
        return None;
    }
    Some(Promoter {
        gene_id: gene.gene_id.clone(),
        symbol: gene.symbol.clone(),
        chrom: gene.chrom.clone(),
        strand: gene.strand,
        start: start0 as u64 + 1,
        end: end0 as u64,
        seq: oriented,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gene(strand: char, tss: u64) -> GeneModel {
        GeneModel {
            gene_id: "G1".into(),
            symbol: "G1".into(),
            chrom: "Chr1".into(),
            start: 1,
            end: 100,
            strand,
            transcript_id: None,
            tss,
        }
    }

    #[test]
    fn plus_strand_takes_upstream() {
        let seq = b"AAAACCCCGGGGTTTT";
        // tss at 1-based 10 -> 0-based 9
        let p = extract_promoter(seq, &gene('+', 10), 4).unwrap();
        assert_eq!(p.start, 6); // 0-based 5..9 -> 1-based 6..9
        assert_eq!(p.end, 9);
        assert_eq!(p.seq.as_slice(), b"CCCG");
    }

    #[test]
    fn minus_strand_reverse_complements_downstream() {
        let seq = b"AAAACCCCGGGGTTTT";
        // tss at 1-based 9 -> 0-based index after TSS is 9; window 9..13
        let p = extract_promoter(seq, &gene('-', 9), 4).unwrap();
        assert_eq!(p.start, 10);
        assert_eq!(p.end, 13);
        // genome 9..13 = GGGT; reverse complement = ACCC
        assert_eq!(p.seq.as_slice(), b"ACCC");
    }
}
