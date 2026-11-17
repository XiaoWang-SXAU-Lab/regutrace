use super::*;
use regutrace_core::{hmmer, tf_classification};

pub(super) fn run(args: &[String]) -> Result<()> {
    let opts = parse_options_with_flags(args, &["--gene-ids"])?;
    for option in opts.keys() {
        if !matches!(
            option.as_str(),
            "--cds" | "--proteins" | "--hmm" | "--hmm-list" | "--rules" | "--gene-map" | "--gene-ids"
                | "--gff" | "--hmmsearch" | "--threads" | "--out"
        ) {
            bail!("unknown option: {option}");
        }
    }
    let required = |name: &str| {
        opts.get(name)
            .map(String::as_str)
            .with_context(|| format!("missing {name}"))
    };
    let (proteins, is_cds) = match (opts.get("--proteins"), opts.get("--cds")) {
        (Some(path), None) => (path.as_str(), false),
        (None, Some(path)) => (path.as_str(), true),
        _ => bail!("supply exactly one of --proteins or --cds"),
    };
    let models = match (opts.get("--hmm"), opts.get("--hmm-list")) {
        (Some(model), None) => vec![PathBuf::from(model)],
        (None, Some(list)) => read_model_list(Path::new(list))?,
        _ => bail!("supply exactly one of --hmm or --hmm-list"),
    };
    let rules = required("--rules")?;
    let out = required("--out")?;
    let gene_map = opts.get("--gene-map");
    let gff_map = opts.get("--gff");
    let gene_ids = opts.get("--gene-ids").is_some_and(|v| {
        matches!(v.to_ascii_lowercase().as_str(), "true" | "1")
    });
    if u8::from(gene_map.is_some()) + u8::from(gff_map.is_some()) + u8::from(gene_ids) != 1 {
        bail!("supply exactly one of --gene-map, --gff or --gene-ids");
    }
    for input in [proteins, rules]
        .into_iter()
        .chain(gene_map.map(String::as_str))
        .chain(gff_map.map(String::as_str))
        .chain(opts.get("--hmm-list").map(String::as_str))
    {
        if same_path(input, out) {
            bail!("input and --out must be different files: {input}");
        }
    }
    for model in &models {
        if same_path(&model.to_string_lossy(), out) {
            bail!("HMM model and --out must be different files: {}", model.display());
        }
    }
    let threads: u32 = opts.get("--threads").map_or("1", String::as_str)
        .parse().context("--threads must be a positive integer")?;
    if threads == 0 {
        bail!("--threads must be a positive integer");
    }
    let pack: tf_classification::RulePack =
        serde_json::from_reader(fasta::open_reader(Path::new(rules))?)
            .with_context(|| format!("read TF rules {rules}"))?;
    let mapping = if let Some(path) = gene_map {
        Some(read_gene_map(Path::new(path))?)
    } else if let Some(path) = gff_map {
        Some(read_gene_map_from_gff(Path::new(path))?)
    } else {
        None
    };
    let out = Path::new(out);
    let tmp = temp_output_path(out);
    fs::create_dir(&tmp).with_context(|| format!("create {}", tmp.display()))?;
    let result = (|| -> Result<usize> {
        let input = tmp.join("proteins.faa");
        let domtblout = tmp.join("domains.txt");
        let ids = if is_cds {
            let translated = tmp.join("translated.faa");
            translate_cds(Path::new(proteins), &translated)?;
            prepare_fasta(&translated, &input, mapping.as_ref())?
        } else {
            prepare_fasta(Path::new(proteins), &input, mapping.as_ref())?
        };
        let program = opts.get("--hmmsearch").map_or("hmmsearch", String::as_str);
        let mut hits = Vec::new();
        for model in &models {
            let search = std::process::Command::new(program)
                .args(["--noali", "--cpu", &threads.to_string(), "--domtblout"])
                .arg(&domtblout)
                .arg(model)
                .arg(&input)
                .stdout(std::process::Stdio::null())
                .output()
                .with_context(|| format!("start {program}; install HMMER or pass --hmmsearch"))?;
            if !search.status.success() {
                bail!(
                    "hmmsearch failed for {} ({}): {}",
                    model.display(), search.status,
                    String::from_utf8_lossy(&search.stderr).trim()
                );
            }
            hits.extend(hmmer::read_hmmsearch_domains(BufReader::new(File::open(&domtblout)?))?);
            fs::remove_file(&domtblout)?;
        }
        if let Some(hit) = hits.iter().find(|hit| !ids.contains(&hit.protein_id)) {
            bail!("hmmsearch returned unknown protein ID: {}", hit.protein_id);
        }
        let (classifications, unknown) = tf_classification::classify(&pack, &hits)?;
        if !unknown.is_empty() {
            eprintln!("models absent from the rule pack: {}", unknown.join(", "));
        }

        let mut genes = BTreeMap::<&str, Vec<_>>::new();
        for classification in &classifications {
            let protein = classification.protein_id;
            let gene = mapping.as_ref().map_or(protein, |map| map[protein].as_str());
            genes.entry(gene).or_default().push(classification);
        }
        let catalogue = tmp.join("tf_genes.tsv");
        write_catalogue(&catalogue, &pack.version, &genes)?;
        publish_temp(&catalogue, out)?;
        Ok(genes.len())
    })();
    for name in ["translated.faa", "proteins.faa", "domains.txt", "tf_genes.tsv"] {
        let _ = fs::remove_file(tmp.join(name));
    }
    let _ = fs::remove_dir(&tmp);
    let count = result?;
    println!("classified {count} TF gene(s) -> {}", out.display());
    Ok(())
}

fn read_model_list(path: &Path) -> Result<Vec<PathBuf>> {
    let base = path.parent().unwrap_or_else(|| Path::new("."));
    let mut models = Vec::new();
    let mut unique = HashSet::new();
    for line in BufReader::new(fasta::open_reader(path)?).lines() {
        let line = line?;
        let line = line.trim_start_matches('\u{feff}').trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let model = base.join(line);
        let canonical = model.canonicalize().with_context(|| format!("open {}", model.display()))?;
        if !unique.insert(canonical) {
            bail!("duplicate HMM model file: {}", model.display());
        }
        models.push(model);
    }
    if models.is_empty() {
        bail!("no HMM model files in {}", path.display());
    }
    Ok(models)
}

fn read_gene_map(path: &Path) -> Result<HashMap<String, String>> {
    let mut reader = ReaderBuilder::new().delimiter(b'\t')
        .from_reader(fasta::open_reader(path)?);
    let header = reader.headers()?;
    let protein = header.iter().position(|name| name == "protein_id")
        .context("gene map requires protein_id")?;
    let gene = header.iter().position(|name| name == "gene_id")
        .context("gene map requires gene_id")?;
    let mut mapping = HashMap::new();
    for row in reader.records() {
        let row = row?;
        let (protein, gene) = (row[protein].trim(), row[gene].trim());
        if protein.is_empty() || gene.is_empty() || gene.contains(char::is_whitespace) {
            bail!("invalid protein_id or gene_id in gene map");
        }
        if mapping.insert(protein.to_owned(), gene.to_owned()).is_some() {
            bail!("duplicate protein_id in gene map: {protein}");
        }
    }
    Ok(mapping)
}

fn read_gene_map_from_gff(path: &Path) -> Result<HashMap<String, String>> {
    parse_gff_gene_map(BufReader::new(fasta::open_reader(path)?).lines())
        .with_context(|| format!("read GFF {}", path.display()))
}

fn parse_gff_gene_map<I>(lines: I) -> Result<HashMap<String, String>>
where
    I: Iterator<Item = std::io::Result<String>>,
{
    let mut gene_ids = HashSet::new();
    let mut transcript_gene = HashMap::<String, String>::new();
    let mut protein_parents = HashMap::<String, Vec<String>>::new();
    let mut protein_locus = HashMap::<String, String>::new();

    for line in lines {
        let line = line?;
        let line = line.trim_end();
        if line.starts_with("##FASTA") {
            break;
        }
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }
        let attrs = parse_gff_attrs(fields[8]);
        match fields[2] {
            "gene" => {
                let id = attrs
                    .get("ID")
                    .or_else(|| attrs.get("Name"))
                    .map(|value| clean_gff_id(value))
                    .filter(|value| !value.is_empty());
                if let Some(id) = id {
                    gene_ids.insert(id);
                }
            }
            "mRNA" | "transcript" => {
                let (Some(id), Some(parents)) = (attrs.get("ID"), attrs.get("Parent")) else {
                    continue;
                };
                let id = clean_gff_id(id);
                for parent in parents.split(',') {
                    let parent = clean_gff_id(parent);
                    if !id.is_empty() && !parent.is_empty() {
                        transcript_gene.insert(id.clone(), parent);
                    }
                }
            }
            "CDS" => {
                let parents: Vec<String> = attrs
                    .get("Parent")
                    .map(|value| {
                        value
                            .split(',')
                            .map(clean_gff_id)
                            .filter(|value| !value.is_empty())
                            .collect()
                    })
                    .unwrap_or_default();
                let locus = attrs.get("locus_tag").map(|value| clean_gff_id(value));
                if let Some(protein_ids) = attrs.get("protein_id") {
                    for protein in protein_ids.split(',') {
                        let protein = protein.trim().trim_matches('"');
                        if protein.is_empty() {
                            continue;
                        }
                        protein_parents
                            .entry(protein.to_owned())
                            .or_default()
                            .extend(parents.iter().cloned());
                        if let Some(locus) = &locus {
                            protein_locus
                                .entry(protein.to_owned())
                                .or_insert_with(|| locus.clone());
                        }
                    }
                } else {
                    // Some annotations (e.g. HF1) omit protein_id and use the
                    // mRNA/transcript ID as the protein FASTA record ID. Treat
                    // each CDS Parent as the protein key and resolve it through
                    // the transcript-to-gene map below.
                    for parent in parents {
                        protein_parents
                            .entry(parent.clone())
                            .or_default()
                            .push(parent.clone());
                        if let Some(locus) = &locus {
                            protein_locus
                                .entry(parent)
                                .or_insert_with(|| locus.clone());
                        }
                    }
                }
            }
            _ => {}
        }
    }

    let proteins: HashSet<&String> = protein_parents
        .keys()
        .chain(protein_locus.keys())
        .collect();
    let mut mapping = HashMap::new();
    for protein in proteins {
        let mut candidates = BTreeSet::new();
        if let Some(parents) = protein_parents.get(protein) {
            for parent in parents {
                if let Some(gene) = transcript_gene.get(parent) {
                    candidates.insert(gene.clone());
                } else if gene_ids.contains(parent) {
                    candidates.insert(parent.clone());
                } else if gene_ids.is_empty() && !parent.is_empty() {
                    // CDS-only GFF: the scan infers genes from the CDS Parent.
                    candidates.insert(parent.clone());
                }
            }
        }
        if candidates.is_empty() {
            if let Some(locus) = protein_locus.get(protein) {
                let prefixed = format!("gene-{locus}");
                if gene_ids.contains(&prefixed) {
                    candidates.insert(prefixed);
                } else if gene_ids.contains(locus) {
                    candidates.insert(locus.clone());
                } else if !gene_ids.is_empty() {
                    candidates.insert(prefixed);
                } else if !locus.is_empty() {
                    candidates.insert(locus.clone());
                }
            }
        }
        if candidates.len() > 1 {
            bail!(
                "protein {protein} maps to multiple genes in GFF: {}",
                candidates.into_iter().collect::<Vec<_>>().join(", ")
            );
        }
        if let Some(gene) = candidates.into_iter().next() {
            mapping.insert(protein.clone(), gene);
        }
    }
    if mapping.is_empty() {
        bail!("no protein-to-gene mappings found in GFF; supply --gene-map or --gene-ids");
    }
    Ok(mapping)
}

fn parse_gff_attrs(attrs: &str) -> HashMap<&str, &str> {
    let mut out = HashMap::new();
    for part in attrs.split(';') {
        let mut kv = part.splitn(2, '=');
        let key = kv.next().map(str::trim).unwrap_or("");
        let value = kv.next().map(str::trim).unwrap_or("");
        if !key.is_empty() && !value.is_empty() {
            out.entry(key).or_insert(value);
        }
    }
    out
}

fn clean_gff_id(value: &str) -> String {
    let value = value.trim().trim_matches('"');
    for prefix in ["gene:", "transcript:", "mRNA:", "rna:"] {
        if let Some(rest) = value.strip_prefix(prefix) {
            return rest.to_string();
        }
    }
    value.to_string()
}

fn prepare_fasta(
    source: &Path,
    destination: &Path,
    mapping: Option<&HashMap<String, String>>,
) -> Result<HashSet<String>> {
    let reader = BufReader::new(fasta::open_reader(source)?);
    let mut writer = BufWriter::new(File::create(destination)?);
    let mut ids = HashSet::new();
    let mut missing = Vec::new();
    let mut needs_sequence = false;
    for line in reader.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            if needs_sequence {
                bail!("protein/CDS FASTA record has no sequence");
            }
            let id = header
                .split_whitespace()
                .next()
                .context("empty protein/CDS FASTA ID")?;
            if !ids.insert(id.to_owned()) {
                bail!("duplicate protein/CDS FASTA ID: {id}");
            }
            if mapping.is_some_and(|map| !map.contains_key(id)) {
                missing.push(id.to_owned());
            }
            needs_sequence = true;
        } else if !line.trim().is_empty() {
            if ids.is_empty() {
                bail!("expected protein/CDS FASTA header before sequence");
            }
            needs_sequence = false;
        }
        writeln!(writer, "{line}")?;
    }
    if ids.is_empty() || needs_sequence {
        bail!("protein/CDS FASTA is empty or ends with an empty record");
    }
    if !missing.is_empty() {
        missing.sort();
        let examples = missing.iter().take(3).cloned().collect::<Vec<_>>().join(", ");
        bail!(
            "protein/CDS FASTA IDs do not match GFF/gene-map IDs ({} unmatched; examples: {}). The FASTA may use transcript/protein IDs instead of gene IDs; ask the AI assistant to generate a protein/CDS-to-gene mapping table, or supply a GFF/gene map with matching protein_id/ID values.",
            missing.len(),
            examples
        );
    }
    writer.flush()?;
    Ok(ids)
}

fn write_catalogue(
    path: &Path,
    version: &str,
    genes: &BTreeMap<&str, Vec<&tf_classification::ProteinClassification<'_>>>,
) -> Result<()> {
    let mut writer = WriterBuilder::new().delimiter(b'\t').from_path(path)?;
    writer.write_record([
        "gene_id", "tf_symbol", "family", "identification_status", "source", "evidence_json",
    ])?;
    for (gene, proteins) in genes {
        let labels: BTreeSet<_> = proteins.iter()
            .flat_map(|protein| protein.families.iter().map(|rule| rule.subfamily.as_str()))
            .collect();
        let status = if labels.len() == 1 { "identified" } else { "ambiguous" };
        let evidence: Vec<_> = proteins.iter().map(|protein| {
            let families: Vec<_> = protein.families.iter().map(|rule| {
                serde_json::json!({"family": rule.family, "subfamily": rule.subfamily})
            }).collect();
            serde_json::json!({
                "protein_id": protein.protein_id,
                "families": families,
                "domains": protein.domains,
            })
        }).collect();
        writer.write_record([
            *gene, *gene, &labels.into_iter().collect::<Vec<_>>().join(";"), status,
            version, &serde_json::to_string(&evidence)?,
        ])?;
    }
    writer.flush()?;
    Ok(())
}

/// Translate complete, oriented CDS with the standard nuclear genetic code.
/// NCBI CDS headers carry the protein accession in `[protein_id=...]`; use it
/// as the record ID so `--gff` protein→gene mapping works without a manual map.
pub(super) fn translate_cds(input: &Path, output: &Path) -> Result<()> {
    let mut reader = BufReader::new(fasta::open_reader(input)?);
    let mut writer = BufWriter::new(File::create(output)?);
    let mut seen = HashSet::new();
    let mut translated_count = 0usize;
    let mut skipped = Vec::new();
    let mut id: Option<String> = None;
    let mut sequence = Vec::new();
    let mut line = String::new();
    loop {
        line.clear();
        let read = reader.read_line(&mut line)?;
        if read == 0 {
            break;
        }
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if let Some(header) = trimmed.strip_prefix('>') {
            if let Some(previous) = id.take() {
                emit_cds_record(&previous, &sequence, &mut writer, &mut seen, &mut translated_count, &mut skipped)?;
            }
            id = Some(cds_record_id(header)?);
            sequence.clear();
        } else if !trimmed.trim().is_empty() {
            if id.is_none() {
                bail!("expected CDS FASTA header before sequence");
            }
            sequence.extend(trimmed.trim().bytes());
        }
    }
    if let Some(previous) = id {
        emit_cds_record(&previous, &sequence, &mut writer, &mut seen, &mut translated_count, &mut skipped)?;
    } else {
        bail!("CDS FASTA is empty");
    }
    writer.flush()?;
    if translated_count == 0 {
        bail!("no translatable CDS records");
    }
    if !skipped.is_empty() {
        eprintln!(
            "skipped {} CDS record(s) with translation errors: {}{}",
            skipped.len(),
            skipped.iter().take(5).cloned().collect::<Vec<_>>().join("; "),
            if skipped.len() > 5 { "; ..." } else { "" }
        );
    }
    Ok(())
}

fn emit_cds_record(
    id: &str,
    sequence: &[u8],
    writer: &mut impl Write,
    seen: &mut HashSet<String>,
    translated_count: &mut usize,
    skipped: &mut Vec<String>,
) -> Result<()> {
    if seen.contains(id) {
        bail!("duplicate CDS ID: {id}");
    }
    match translate_cds_record(id, sequence, writer, seen) {
        Ok(()) => {
            *translated_count += 1;
            Ok(())
        }
        Err(error) => {
            skipped.push(format!("{id}: {error}"));
            Ok(())
        }
    }
}

fn cds_record_id(header: &str) -> Result<String> {
    let first = header.split_whitespace().next().unwrap_or("").trim();
    if first.is_empty() {
        bail!("empty CDS FASTA ID");
    }
    if let Some(position) = header.find("protein_id=") {
        let rest = &header[position + "protein_id=".len()..];
        let end = rest
            .find(|c: char| c.is_whitespace() || c == ']' || c == ';')
            .unwrap_or(rest.len());
        let protein = rest[..end].trim();
        if !protein.is_empty() {
            return Ok(protein.to_owned());
        }
    }
    Ok(first.to_owned())
}

fn translate_cds_record(
    id: &str,
    sequence: &[u8],
    writer: &mut impl Write,
    seen: &mut HashSet<String>,
) -> Result<()> {
    if !seen.insert(id.to_owned()) {
        bail!("duplicate CDS ID: {id}");
    }
    if sequence.is_empty() || sequence.len() % 3 != 0 {
        bail!("CDS {id} must have a positive length divisible by three");
    }
    let mut translated = Vec::with_capacity(sequence.len() / 3);
    for codon in sequence.chunks_exact(3) {
        let base = |b: u8| b"TCAG".iter().position(|&x| x == b);
        let amino = match (base(codon[0]), base(codon[1]), base(codon[2])) {
            (Some(a), Some(b), Some(c)) => {
                b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[a * 16 + b * 4 + c]
            }
            _ if codon.iter().all(|b| b"TCAGRYSWKMBDHVN".contains(b)) => b'X',
            _ => bail!("invalid nucleotide in CDS {id}"),
        };
        translated.push(amino);
    }
    if translated.last() == Some(&b'*') {
        translated.pop();
    }
    if translated.is_empty() || translated.contains(&b'*') {
        bail!("CDS {id} has an internal stop or no translated residues");
    }
    writeln!(writer, ">{id}")?;
    writer.write_all(&translated)?;
    writeln!(writer)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{cds_record_id, parse_gff_gene_map, prepare_fasta};

    fn parse(lines: &[&str]) -> anyhow::Result<std::collections::HashMap<String, String>> {
        parse_gff_gene_map(lines.iter().map(|line| Ok((*line).to_owned())))
    }

    #[test]
    fn prepare_fasta_reports_unmatched_ids_with_ai_guidance() {
        let dir = std::env::temp_dir();
        let source = dir.join(format!("regutrace_prepare_{}.faa", std::process::id()));
        let destination = dir.join(format!("regutrace_prepare_out_{}.faa", std::process::id()));
        std::fs::write(&source, ">P1\nACDE\n>P2\nFGHI\n").unwrap();
        let mapping = std::collections::HashMap::new();
        let err = prepare_fasta(&source, &destination, Some(&mapping))
            .unwrap_err()
            .to_string();
        assert!(err.contains("do not match"), "{err}");
        assert!(err.contains("AI assistant"), "{err}");
        assert!(err.contains("P1"), "{err}");
        let _ = std::fs::remove_file(source);
        let _ = std::fs::remove_file(destination);
    }

    #[test]
    fn cds_uses_protein_id_when_present() {
        assert_eq!(
            cds_record_id("lcl|NC_1_cds_NP_009332.1_1 [protein_id=NP_009332.1] [locus_tag=YAL068C]").unwrap(),
            "NP_009332.1"
        );
        assert_eq!(cds_record_id("AT1G01010.1").unwrap(), "AT1G01010.1");
        assert_eq!(
            cds_record_id("tx1 [protein_id=P1;locus_tag=G1]").unwrap(),
            "P1"
        );
    }

    #[test]
    fn resolves_transcript_parent_to_gene_feature() {
        let mapping = parse(&[
            "##gff-version 3",
            "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene-YAL001W;Name=YAL001W",
            "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=rna-XM_1;Parent=gene-YAL001W",
            "chr1\tsrc\tCDS\t1\t100\t.\t+\t0\tID=cds-NP_1;Parent=rna-XM_1;protein_id=NP_1;locus_tag=YAL001W",
        ]).unwrap();
        assert_eq!(mapping.get("NP_1").map(String::as_str), Some("gene-YAL001W"));
    }

    #[test]
    fn resolves_transcript_parent_without_protein_id() {
        let mapping = parse(&[
            "##gff-version 3",
            "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene-HF1",
            "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=HF1.1;Parent=gene-HF1",
            "chr1\tsrc\tCDS\t1\t100\t.\t+\t0\tParent=HF1.1",
        ]).unwrap();
        assert_eq!(mapping.get("HF1.1").map(String::as_str), Some("gene-HF1"));
    }

    #[test]
    fn resolves_cds_parent_that_is_a_gene() {
        let mapping = parse(&[
            "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=G1",
            "chr1\tsrc\tCDS\t1\t100\t.\t+\t0\tParent=G1;protein_id=P1",
        ]).unwrap();
        assert_eq!(mapping.get("P1").map(String::as_str), Some("G1"));
    }

    #[test]
    fn keeps_cds_parent_for_cds_only_gff() {
        let mapping = parse(&[
            "chr1\tsrc\tCDS\t1\t100\t.\t+\t0\tParent=tx1;protein_id=P1",
        ]).unwrap();
        assert_eq!(mapping.get("P1").map(String::as_str), Some("tx1"));
    }

    #[test]
    fn falls_back_to_prefixed_locus_tag() {
        let mapping = parse(&[
            "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene-YAL003W",
            "chr1\tsrc\tCDS\t1\t100\t.\t+\t0\tParent=tx3;protein_id=P3;locus_tag=YAL003W",
        ]).unwrap();
        assert_eq!(mapping.get("P3").map(String::as_str), Some("gene-YAL003W"));
    }

    #[test]
    fn rejects_ambiguous_protein_mapping() {
        let err = parse(&[
            "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene-A",
            "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene-B",
            "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=tx1;Parent=gene-A",
            "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=tx2;Parent=gene-B",
            "chr1\tsrc\tCDS\t1\t100\t.\t+\t0\tParent=tx1,tx2;protein_id=P1",
        ]).unwrap_err();
        assert!(err.to_string().contains("multiple genes"));
    }
}
