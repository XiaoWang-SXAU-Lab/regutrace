//! Minimal FASTA reader (plain or gzip) with chromosome-name fallback.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;

pub struct Fasta {
    pub seqs: HashMap<String, Vec<u8>>,
    pub order: Vec<String>,
}

pub fn open_reader(path: &Path) -> Result<Box<dyn Read>> {
    let mut file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut magic = [0u8; 2];
    let n = file.read(&mut magic)?;
    file.seek(SeekFrom::Start(0))?;
    let gz_by_magic = n == 2 && magic == [0x1f, 0x8b];
    let gz_by_extension = path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);
    if gz_by_magic || gz_by_extension {
        Ok(Box::new(MultiGzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

/// Split a sequence name into a normalized prefix family and numeric core.
/// Examples: chr01 / chromosome1 / CHR1 -> ("chr", "1");
/// scaffold_001 -> ("scaffold", "1"); 1 -> (None, "1").
fn split_chrom(name: &str) -> (Option<&'static str>, String) {
    let lower = name.trim().to_ascii_lowercase();
    for (prefixes, family) in [
        (&["chromosome", "chrom", "chr"][..], "chr"),
        (&["scaffold", "scaf"][..], "scaffold"),
        (&["contig", "ctg"][..], "contig"),
        (&["linkage_group", "linkagegroup", "lg"][..], "lg"),
    ] {
        for prefix in prefixes {
            if let Some(rest) = lower.strip_prefix(prefix) {
                let rest = rest.trim_start_matches(|c: char| c == '_' || c == '-' || c == ' ');
                let core = rest.trim_start_matches('0');
                return (Some(family), if core.is_empty() { "0" } else { core }.to_string());
            }
        }
    }
    let core = lower.trim_start_matches('0');
    (None, if core.is_empty() { "0" } else { core }.to_string())
}

/// Score alias candidates. Lower is better; None means the prefix families are
/// different (for example scaffold_1 must not silently match chr1).
fn chrom_match_score(query: Option<&'static str>, key: Option<&'static str>) -> Option<u8> {
    match (query, key) {
        (Some(a), Some(b)) if a == b => Some(0),
        (None, None) => Some(1),
        (Some("chr"), None) | (None, Some("chr")) => Some(2),
        (Some(_), None) | (None, Some(_)) => Some(3),
        _ => None,
    }
}

impl Fasta {
    /// Resolve common chromosome naming differences: case, `chr`/`chromosome`,
    /// `scaffold`/`scaf`, `contig`/`ctg`, and zero padding (`chr1` vs `chr01`).
    pub fn seq(&self, chrom: &str) -> Option<&Vec<u8>> {
        if let Some(s) = self.seqs.get(chrom) {
            return Some(s);
        }
        if let Some((_, s)) = self
            .seqs
            .iter()
            .find(|(name, _)| name.eq_ignore_ascii_case(chrom))
        {
            return Some(s);
        }
        let (query_family, query_core) = split_chrom(chrom);
        let mut keys: Vec<&String> = self.seqs.keys().collect();
        keys.sort();
        let mut best: Option<(u8, &Vec<u8>)> = None;
        for key in keys {
            let (key_family, key_core) = split_chrom(key);
            if key_core != query_core {
                continue;
            }
            let Some(score) = chrom_match_score(query_family, key_family) else {
                continue;
            };
            if best.as_ref().is_none_or(|(best_score, _)| score < *best_score) {
                best = Some((score, self.seqs.get(key).unwrap()));
            }
        }
        best.map(|(_, seq)| seq)
    }
}

pub fn read_fasta(path: &Path) -> Result<Fasta> {
    let mut reader = BufReader::new(open_reader(path)?);
    let mut seqs = HashMap::new();
    let mut order = Vec::new();
    let mut name: Option<String> = None;
    let mut buf: Vec<u8> = Vec::new();
    let mut line = String::new();

    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        let trimmed = line.trim_end_matches(|c| c == '\r' || c == '\n');
        if let Some(rest) = trimmed.strip_prefix('>') {
            if let Some(prev) = name.take() {
                if buf.is_empty() {
                    bail!("FASTA record {prev} has no sequence");
                }
                if seqs.contains_key(&prev) {
                    bail!("duplicate FASTA record: {prev}");
                }
                order.push(prev.clone());
                seqs.insert(prev, std::mem::take(&mut buf));
            }
            let id = rest.split_whitespace().next().unwrap_or("").to_string();
            if id.is_empty() {
                bail!("FASTA header without an identifier");
            }
            name = Some(id);
        } else if name.is_some() {
            for b in trimmed.bytes() {
                if b.is_ascii_whitespace() {
                    continue;
                }
                buf.push(b.to_ascii_uppercase());
            }
        }
    }

    if let Some(prev) = name {
        if buf.is_empty() {
            bail!("FASTA record {prev} has no sequence");
        }
        if seqs.contains_key(&prev) {
            bail!("duplicate FASTA record: {prev}");
        }
        order.push(prev.clone());
        seqs.insert(prev, buf);
    }
    if seqs.is_empty() {
        bail!("no FASTA records in {}", path.display());
    }
    Ok(Fasta { seqs, order })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn reads_multiline_fasta() {
        let path = std::env::temp_dir().join(format!("regutrace_fasta_{}.fa", std::process::id()));
        {
            let mut f = File::create(&path).unwrap();
            writeln!(f, ">Chr1 test").unwrap();
            writeln!(f, "ACGT").unwrap();
            writeln!(f, "acgt").unwrap();
        }
        let fa = read_fasta(&path).unwrap();
        assert_eq!(fa.seq("Chr1").unwrap().as_slice(), b"ACGTACGT");
        assert_eq!(fa.seq("1").unwrap().as_slice(), b"ACGTACGT");
        assert_eq!(fa.seq("chr1").unwrap().as_slice(), b"ACGTACGT");
        assert_eq!(fa.seq("CHR1").unwrap().as_slice(), b"ACGTACGT");
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn resolves_chromosome_aliases() {
        let mut seqs = HashMap::new();
        seqs.insert("chr1".to_string(), b"AAA".to_vec());
        seqs.insert("scaffold_1".to_string(), b"CCC".to_vec());
        let fa = Fasta {
            seqs,
            order: vec![],
        };
        assert_eq!(fa.seq("Chr1").unwrap().as_slice(), b"AAA");
        assert_eq!(fa.seq("CHR01").unwrap().as_slice(), b"AAA");
        assert_eq!(fa.seq("chromosome1").unwrap().as_slice(), b"AAA");
        assert_eq!(fa.seq("1").unwrap().as_slice(), b"AAA");
        assert_eq!(fa.seq("scaffold1").unwrap().as_slice(), b"CCC");
        assert!(fa.seq("contig1").is_none());
    }

    #[test]
    fn reads_gzip_without_gz_extension() {
        use flate2::Compression;
        use flate2::write::GzEncoder;
        let path = std::env::temp_dir().join(format!(
            "regutrace_fasta_gz_{}.fa",
            std::process::id()
        ));
        {
            let f = File::create(&path).unwrap();
            let mut enc = GzEncoder::new(f, Compression::default());
            enc.write_all(b">Chr1\nACGT\n").unwrap();
            enc.finish().unwrap();
        }
        let fa = read_fasta(&path).unwrap();
        assert_eq!(fa.seq("Chr1").unwrap().as_slice(), b"ACGT");
        let _ = std::fs::remove_file(path);
    }
}
