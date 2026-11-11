//! Minimal GFF3 reader for gene/transcript models.
//!
//! Coordinates are kept 1-based inclusive as in GFF3. `tss` is the 1-based
//! coordinate of the transcription start site on the gene's strand.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::fasta::open_reader;

#[derive(Clone, Debug)]
pub struct GeneModel {
    pub gene_id: String,
    pub symbol: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub transcript_id: Option<String>,
    pub tss: u64,
}

#[derive(Clone, Debug)]
struct GeneInfo {
    symbol: String,
    chrom: String,
    start: u64,
    end: u64,
    strand: char,
}

#[derive(Clone, Debug)]
struct TxInfo {
    transcript_id: String,
    start: u64,
    end: u64,
}

fn attr(attrs: &str, key: &str) -> Option<String> {
    for part in attrs.split(';') {
        let mut kv = part.splitn(2, '=');
        let k = kv.next()?.trim();
        if k != key {
            continue;
        }
        let v = kv.next()?.trim().trim_matches('"');
        if !v.is_empty() {
            return Some(v.to_string());
        }
    }
    None
}

fn clean_id(value: &str) -> String {
    let v = value.trim();
    for prefix in ["gene:", "transcript:", "mRNA:", "rna:"] {
        if let Some(rest) = v.strip_prefix(prefix) {
            return rest.to_string();
        }
    }
    v.to_string()
}

pub fn read_gff(path: &Path) -> Result<Vec<GeneModel>> {
    let reader = BufReader::new(open_reader(path)?);
    let mut genes: HashMap<String, GeneInfo> = HashMap::new();
    let mut tx: HashMap<String, Vec<TxInfo>> = HashMap::new();

    for (lineno, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("read GFF line {}", lineno + 1))?;
        let line = line.trim_end_matches(|c| c == '\r' || c == '\n');
        if line.is_empty() || line.starts_with('#') {
            if line.starts_with("##FASTA") {
                break;
            }
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            bail!("GFF line {} has {} fields", lineno + 1, fields.len());
        }
        let chrom = fields[0].to_string();
        let kind = fields[2];
        let start: u64 = fields[3]
            .parse()
            .with_context(|| format!("GFF start at line {}", lineno + 1))?;
        let end: u64 = fields[4]
            .parse()
            .with_context(|| format!("GFF end at line {}", lineno + 1))?;
        let strand = fields[6].chars().next().unwrap_or('+');
        let attrs = fields[8];

        if kind == "gene" {
            let Some(gene_id_raw) = attr(attrs, "ID").or_else(|| attr(attrs, "Name")) else {
                continue;
            };
            let gene_id = clean_id(&gene_id_raw);
            let symbol = attr(attrs, "Name")
                .or_else(|| attr(attrs, "gene_name"))
                .unwrap_or_else(|| gene_id.clone());
            genes.insert(
                gene_id.clone(),
                GeneInfo {
                    symbol,
                    chrom,
                    start,
                    end,
                    strand,
                },
            );
        } else if kind == "mRNA" || kind == "transcript" {
            let Some(parent_raw) = attr(attrs, "Parent") else {
                continue;
            };
            let parent = clean_id(&parent_raw);
            let transcript_id = attr(attrs, "ID")
                .map(|v| clean_id(&v))
                .unwrap_or_else(|| format!("{parent}.t{}", tx.get(&parent).map_or(1, |v| v.len() + 1)));
            tx.entry(parent).or_default().push(TxInfo {
                transcript_id,
                start,
                end,
            });
        }
    }

    // GFF3 without explicit gene features: use transcript parent as the gene.
    for (parent, txs) in &tx {
        if genes.contains_key(parent) {
            continue;
        }
        if let Some(first) = txs.first() {
            genes.insert(
                parent.clone(),
                GeneInfo {
                    symbol: parent.clone(),
                    chrom: String::new(),
                    start: first.start,
                    end: first.end,
                    strand: '+',
                },
            );
        }
    }

    let mut out = Vec::with_capacity(genes.len());
    for (gene_id, info) in genes {
        let best = tx
            .get(&gene_id)
            .and_then(|v| v.iter().max_by_key(|t| t.end.saturating_sub(t.start)));
        let (transcript_id, tss) = match best {
            Some(t) => (
                Some(t.transcript_id.clone()),
                if info.strand == '-' { t.end } else { t.start },
            ),
            None => (
                None,
                if info.strand == '-' { info.end } else { info.start },
            ),
        };
        out.push(GeneModel {
            gene_id,
            symbol: info.symbol,
            chrom: info.chrom,
            start: info.start,
            end: info.end,
            strand: info.strand,
            transcript_id,
            tss,
        });
    }
    out.sort_by(|a, b| a.gene_id.cmp(&b.gene_id));
    if out.is_empty() {
        bail!("no gene features found in {}", path.display());
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn parses_plus_and_minus_genes() {
        let path = std::env::temp_dir().join(format!("regutrace_gff_{}.gff3", std::process::id()));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "##gff-version 3").unwrap();
            writeln!(f, "Chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene:AT1G00010;Name=PLUS").unwrap();
            writeln!(f, "Chr1\ttest\tmRNA\t100\t500\t.\t+\t.\tID=transcript:AT1G00010.1;Parent=gene:AT1G00010").unwrap();
            writeln!(f, "Chr2\ttest\tgene\t900\t1200\t.\t-\t.\tID=gene:AT2G00020;Name=MINUS").unwrap();
            writeln!(f, "Chr2\ttest\tmRNA\t900\t1100\t.\t-\t.\tID=transcript:AT2G00020.1;Parent=gene:AT2G00020").unwrap();
        }
        let genes = read_gff(&path).unwrap();
        let plus = genes.iter().find(|g| g.gene_id == "AT1G00010").unwrap();
        assert_eq!(plus.strand, '+');
        assert_eq!(plus.tss, 100);
        let minus = genes.iter().find(|g| g.gene_id == "AT2G00020").unwrap();
        assert_eq!(minus.strand, '-');
        assert_eq!(minus.tss, 1100);
        let _ = std::fs::remove_file(path);
    }
}
