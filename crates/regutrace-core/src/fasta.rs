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

impl Fasta {
    /// Resolve common plant chromosome naming differences (`1`, `Chr1`,
    /// `chr1`, `CHR1`).
    pub fn seq(&self, chrom: &str) -> Option<&Vec<u8>> {
        if let Some(s) = self.seqs.get(chrom) {
            return Some(s);
        }
        let core = chrom
            .strip_prefix("Chr")
            .or_else(|| chrom.strip_prefix("chr"))
            .or_else(|| chrom.strip_prefix("CHR"))
            .unwrap_or(chrom);
        for c in [
            format!("Chr{core}"),
            format!("chr{core}"),
            format!("CHR{core}"),
            core.to_string(),
        ] {
            if c != chrom {
                if let Some(s) = self.seqs.get(&c) {
                    return Some(s);
                }
            }
        }
        None
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
