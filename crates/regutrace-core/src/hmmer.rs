use std::io::BufRead;

use anyhow::{ensure, Context, Result};

/// One reported domain from `hmmsearch --domtblout`.
#[derive(Debug, PartialEq, serde::Serialize)]
pub struct DomainHit {
    pub protein_id: String,
    pub model_name: String,
    pub model_accession: Option<String>,
    pub model_length: u32,
    pub full_evalue: f64,
    pub domain_score: f64,
    pub conditional_evalue: f64,
    pub independent_evalue: f64,
    pub model_start: u32,
    pub model_end: u32,
    pub protein_start: u32,
    pub protein_end: u32,
}

/// Preserve repeated domains and their order for downstream family rules.
pub fn read_hmmsearch_domains(reader: impl BufRead) -> Result<Vec<DomainHit>> {
    let mut hits = Vec::new();
    for (index, line) in reader.lines().enumerate() {
        let context = || format!("hmmsearch domtblout line {}", index + 1);
        let line = line.with_context(context)?;
        let row = line.trim();
        if row.is_empty() || row.starts_with('#') {
            continue;
        }
        hits.push(parse_domain(row).with_context(context)?);
    }
    Ok(hits)
}

fn parse_domain(row: &str) -> Result<DomainHit> {
    let fields: Vec<_> = row.split_whitespace().take(22).collect();
    ensure!(fields.len() == 22, "expected 22 fixed fields");

    let model_length = fields[5].parse()?;
    let protein_length: u32 = fields[2].parse()?;
    let hit = DomainHit {
        protein_id: fields[0].to_owned(),
        model_name: fields[3].to_owned(),
        model_accession: (fields[4] != "-").then(|| fields[4].to_owned()),
        model_length,
        full_evalue: parse_evalue(fields[6])?,
        domain_score: fields[13].parse()?,
        conditional_evalue: parse_evalue(fields[11])?,
        independent_evalue: parse_evalue(fields[12])?,
        model_start: fields[15].parse()?,
        model_end: fields[16].parse()?,
        protein_start: fields[17].parse()?,
        protein_end: fields[18].parse()?,
    };
    ensure!(hit.domain_score.is_finite(), "non-finite domain score");
    ensure!(
        hit.model_start > 0 && hit.model_start <= hit.model_end && hit.model_end <= model_length,
        "invalid model coordinates"
    );
    ensure!(
        hit.protein_start > 0
            && hit.protein_start <= hit.protein_end
            && hit.protein_end <= protein_length,
        "invalid protein coordinates"
    );
    Ok(hit)
}

fn parse_evalue(field: &str) -> Result<f64> {
    let value: f64 = field.parse()?;
    ensure!(value.is_finite() && value >= 0.0, "invalid E-value");
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn preserves_repeated_domains_and_distinct_evalues() {
        let input = "# hmmsearch output\n\
            TX1 - 300 AP2 PF00847.21 60 1e-19 65.0 0.1 1 2 1e-9 1e-8 38.0 0.0 1 55 10 64 8 65 0.93 protein description\n\
            TX1 - 300 AP2 PF00847.21 60 1e-19 65.0 0.1 2 2 1e-7 1e-6 31.0 0.0 2 57 110 165 108 165 0.91 protein description\n";
        let hits = read_hmmsearch_domains(input.as_bytes()).unwrap();
        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0].protein_id, "TX1");
        assert_eq!(hits[0].model_accession.as_deref(), Some("PF00847.21"));
        assert_eq!(hits[0].full_evalue, 1e-19);
        assert_eq!(hits[0].conditional_evalue, 1e-9);
        assert_eq!(hits[0].independent_evalue, 1e-8);
        assert_eq!(hits[1].protein_start, 110);
    }

    #[test]
    fn reports_the_line_of_a_truncated_row() {
        let err = read_hmmsearch_domains(b"# header\nTX1 - 300\n".as_slice()).unwrap_err();
        assert_eq!(err.to_string(), "hmmsearch domtblout line 2");
    }
}
