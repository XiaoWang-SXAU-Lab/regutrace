//! End-to-end check for the regutrace binary on a synthetic genome.

use std::path::PathBuf;
use std::process::Command;

fn write(path: &PathBuf, text: &str) {
    std::fs::write(path, text).unwrap();
}

#[test]
fn cli_end_to_end() {
    let dir = std::env::temp_dir().join(format!("regutrace_cli_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let motif = "ACGTAC";
    let mut genome = vec![b'A'; 3000];
    genome[150..156].copy_from_slice(motif.as_bytes());
    // Gene 2 on the minus strand: put the reverse complement in the promoter.
    let rc: Vec<u8> = motif
        .bytes()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect();
    genome[800..806].copy_from_slice(&rc);

    let genome_path = dir.join("genome.fa");
    write(&genome_path, &format!(">Chr1\n{}\n", String::from_utf8(genome).unwrap()));

    let gff_path = dir.join("annotation.gff3");
    write(
        &gff_path,
        "##gff-version 3\n\
Chr1\ttest\tgene\t201\t400\t.\t+\t.\tID=gene:G1;Name=GENE1\n\
Chr1\ttest\tmRNA\t201\t400\t.\t+\t.\tID=transcript:G1.t1;Parent=gene:G1\n\
Chr1\ttest\tgene\t601\t800\t.\t-\t.\tID=gene:G2;Name=GENE2\n\
Chr1\ttest\tmRNA\t601\t800\t.\t-\t.\tID=transcript:G2.t1;Parent=gene:G2\n",
    );

    let motif_path = dir.join("motifs.txt");
    let mut pwm = String::from(">MA0001.1 TEST\n");
    for base in ["A", "C", "G", "T"] {
        let values: Vec<String> = motif
            .chars()
            .map(|m| if m.to_string() == base { "10".to_string() } else { "0".to_string() })
            .collect();
        pwm.push_str(&format!("{base} [ {} ]\n", values.join(" ")));
    }
    write(&motif_path, &pwm);

    let map_path = dir.join("tf_map.tsv");
    write(
        &map_path,
        "tf_id\ttf_symbol\tfamily\tmotif_id\nAT1G00010\tTF1\tHSF\tMA0001.1\n",
    );
    let anno_path = dir.join("annotation.tsv");
    write(&anno_path, "gene_id\tsymbol\nG1\tGENE1\nG2\tGENE2\n");
    let user_path = dir.join("user.tsv");
    write(
        &user_path,
        "tf\ttarget\tevidence_type\tevidence_level\tmethod\tsource\tpmid\n\
AT1G00010\tG2\tpoint_validation_emsa\tA\tEMSA\tPMID:123\t123\n",
    );
    let out_path = dir.join("out.tsv");
    let scan_report = dir.join("scan-report.json");
    let progress_path = dir.join("progress.json");

    let bin = env!("CARGO_BIN_EXE_regutrace");
    let output = Command::new(bin)
        .args([
            "scan",
            "--genome",
            genome_path.to_str().unwrap(),
            "--gff",
            gff_path.to_str().unwrap(),
            "--motifs",
            motif_path.to_str().unwrap(),
            "--tf-map",
            map_path.to_str().unwrap(),
            "--promoter-len",
            "100",
            "--downstream",
            "100",
            "--pvalue",
            "1e-3",
            "--annotation",
            anno_path.to_str().unwrap(),
            "--user-list",
            user_path.to_str().unwrap(),
            "--out",
            out_path.to_str().unwrap(),
            "--report-json",
            scan_report.to_str().unwrap(),
            "--progress-json",
            progress_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(&out_path).unwrap();
    assert!(text.contains("prediction_motif"), "{text}");
    assert!(text.contains("point_validation_emsa"), "{text}");
    assert!(text.contains("GENE1"), "{text}");
    assert!(text.contains("GENE2"), "{text}");
    assert!(text.lines().count() >= 3, "{text}");

    let scan_json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&scan_report).unwrap()).unwrap();
    assert_eq!(scan_json["command"], "scan");
    assert_eq!(scan_json["summary"]["errors"], 0);
    assert!(scan_json["diagnostics"]
        .as_array()
        .unwrap()
        .iter()
        .any(|d| d["code"] == "SCAN_SUMMARY"));

    let progress_json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&progress_path).unwrap()).unwrap();
    assert_eq!(progress_json["status"], "done");
    assert_eq!(progress_json["stage"], "done");
    assert_eq!(progress_json["progress"], 1.0);
    assert_eq!(progress_json["done"], 1);
    assert_eq!(progress_json["total"], 1);

    let rows: Vec<Vec<&str>> = text.lines().skip(1).map(|l| l.split('\t').collect()).collect();
    // Columns: 0 tf_id, 1 tf_symbol, 2 target_id, 3 target_symbol, 4 evidence_type,
    // 5 evidence_level, 10 promoter_start, 11 promoter_end, 13 hit_start,
    // 14 hit_end, 15 strand.
    let g1 = rows
        .iter()
        .find(|r| r[2] == "G1" && r[4] == "prediction_motif")
        .expect("G1 prediction");
    assert_eq!(g1[10], "101");
    assert_eq!(g1[11], "300");
    assert_eq!(g1[13], "51");
    assert_eq!(g1[14], "56");
    assert_eq!(g1[15], "+");
    let g2 = rows
        .iter()
        .find(|r| r[2] == "G2" && r[4] == "prediction_motif")
        .expect("G2 prediction");
    assert_eq!(g2[10], "701");
    assert_eq!(g2[11], "900");
    assert_eq!(g2[13], "95");
    assert_eq!(g2[14], "100");
    assert_eq!(g2[15], "+");

    // Without --tf-map, TF IDs are derived from the motif header name.
    let auto_out = dir.join("auto.tsv");
    let auto = Command::new(bin)
        .args([
            "scan",
            "--genome",
            genome_path.to_str().unwrap(),
            "--gff",
            gff_path.to_str().unwrap(),
            "--motifs",
            motif_path.to_str().unwrap(),
            "--promoter-len",
            "100",
            "--downstream",
            "100",
            "--pvalue",
            "1e-3",
            "--out",
            auto_out.to_str().unwrap(),
            "--threads",
            "1",
        ])
        .output()
        .unwrap();
    assert!(
        auto.status.success(),
        "auto stderr: {}",
        String::from_utf8_lossy(&auto.stderr)
    );
    let auto_text = std::fs::read_to_string(&auto_out).unwrap();
    assert!(
        auto_text.lines().any(|l| l.starts_with("TEST\tTEST\tG1\t")),
        "{auto_text}"
    );

    // Streaming output is used by the dashboard for whole-genome scans.
    let stream_out = dir.join("stream.tsv");
    let stream_report = dir.join("stream-report.json");
    let stream = Command::new(bin)
        .args([
            "scan",
            "--genome",
            genome_path.to_str().unwrap(),
            "--gff",
            gff_path.to_str().unwrap(),
            "--motifs",
            motif_path.to_str().unwrap(),
            "--tf-map",
            map_path.to_str().unwrap(),
            "--promoter-len",
            "100",
            "--downstream",
            "100",
            "--pvalue",
            "1e-3",
            "--out",
            stream_out.to_str().unwrap(),
            "--report-json",
            stream_report.to_str().unwrap(),
            "--stream",
            "--threads",
            "2",
        ])
        .output()
        .unwrap();
    assert!(
        stream.status.success(),
        "stream stderr: {}",
        String::from_utf8_lossy(&stream.stderr)
    );
    let stream_text = std::fs::read_to_string(&stream_out).unwrap();
    assert!(stream_text.contains("prediction_motif"), "{stream_text}");
    assert!(stream_text.lines().count() >= 3, "{stream_text}");
    let stream_json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&stream_report).unwrap()).unwrap();
    let summary = stream_json["diagnostics"]
        .as_array()
        .unwrap()
        .iter()
        .find(|d| d["code"] == "SCAN_SUMMARY")
        .expect("SCAN_SUMMARY");
    assert_eq!(summary["context"]["motif_hits"], "2");
    assert_eq!(summary["context"]["total_edges"], "2");

    let validate_report = dir.join("validate-report.json");
    let validate = Command::new(bin)
        .args([
            "validate",
            "--genome",
            genome_path.to_str().unwrap(),
            "--gff",
            gff_path.to_str().unwrap(),
            "--motifs",
            motif_path.to_str().unwrap(),
            "--tf-map",
            map_path.to_str().unwrap(),
            "--annotation",
            anno_path.to_str().unwrap(),
            "--user-list",
            user_path.to_str().unwrap(),
            "--promoter-len",
            "100",
            "--downstream",
            "100",
            "--report-json",
            validate_report.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        validate.status.success(),
        "validate stderr: {}",
        String::from_utf8_lossy(&validate.stderr)
    );
    let validate_json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&validate_report).unwrap()).unwrap();
    assert_eq!(validate_json["command"], "validate");
    assert_eq!(validate_json["summary"]["errors"], 0);
    assert!(validate_json["diagnostics"]
        .as_array()
        .unwrap()
        .iter()
        .any(|d| d["severity"] == "warning"));

    let strict = Command::new(bin)
        .args([
            "validate",
            "--genome",
            genome_path.to_str().unwrap(),
            "--gff",
            gff_path.to_str().unwrap(),
            "--motifs",
            motif_path.to_str().unwrap(),
            "--promoter-len",
            "100",
            "--strict",
        ])
        .output()
        .unwrap();
    assert!(!strict.status.success());
    assert!(
        String::from_utf8_lossy(&strict.stderr).contains("STRICT_WARNINGS"),
        "{}",
        String::from_utf8_lossy(&strict.stderr)
    );

    let arg_report = dir.join("argfail.json");
    let argfail = Command::new(bin)
        .args([
            "scan",
            "--pvalue",
            "2",
            "--report-json",
            arg_report.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!argfail.status.success());
    let arg_json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&arg_report).unwrap()).unwrap();
    assert_eq!(arg_json["command"], "scan");
    assert!(arg_json["diagnostics"]
        .as_array()
        .unwrap()
        .iter()
        .any(|d| d["code"] == "SCAN_FAILED"));

    let db_path = dir.join("evidence.db");
    let import = Command::new(bin)
        .args([
            "import",
            "--db",
            db_path.to_str().unwrap(),
            "--evidence",
            out_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        import.status.success(),
        "import stderr: {}",
        String::from_utf8_lossy(&import.stderr)
    );

    let stats = Command::new(bin)
        .args(["stats", "--db", db_path.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(stats.status.success());
    let stats_text = String::from_utf8_lossy(&stats.stdout);
    assert!(stats_text.contains("evidence_level"), "{stats_text}");
    assert!(stats_text.contains("prediction_motif"), "{stats_text}");

    let export_path = dir.join("export.tsv");
    let export = Command::new(bin)
        .args([
            "export",
            "--db",
            db_path.to_str().unwrap(),
            "--out",
            export_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        export.status.success(),
        "export stderr: {}",
        String::from_utf8_lossy(&export.stderr)
    );
    let exported = std::fs::read_to_string(&export_path).unwrap();
    assert!(exported.contains("prediction_motif"), "{exported}");
    assert!(exported.contains("point_validation_emsa"), "{exported}");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn cds_only_gff_without_biotype_respects_pack_filter_policy() {
    let dir = std::env::temp_dir().join(format!(
        "regutrace_cli_biotype_{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let genome_path = dir.join("genome.fa");
    write(&genome_path, &format!(">Chr1\n{}\n", "A".repeat(2000)));

    // CDS-only annotation: no gene/mRNA features and no gene_biotype attr.
    let gff_path = dir.join("cds-only.gff3");
    write(
        &gff_path,
        "##gff-version 3\n\
Chr1\ttest\tCDS\t201\t400\t.\t+\t0\tParent=geneA\n\
Chr1\ttest\tCDS\t601\t800\t.\t+\t0\tParent=geneB\n",
    );

    let motif_path = dir.join("motifs.txt");
    write(
        &motif_path,
        ">MA_TEST.1 ALL_A\n\
A [ 10 10 10 10 ]\n\
C [ 0 0 0 0 ]\n\
G [ 0 0 0 0 ]\n\
T [ 0 0 0 0 ]\n",
    );

    let bin = env!("CARGO_BIN_EXE_regutrace");
    let explicit_out = dir.join("explicit.tsv");
    let explicit = Command::new(bin)
        .args([
            "scan",
            "--genome",
            genome_path.to_str().unwrap(),
            "--gff",
            gff_path.to_str().unwrap(),
            "--motifs",
            motif_path.to_str().unwrap(),
            "--gene-biotype",
            "protein_coding",
            "--promoter-len",
            "100",
            "--pvalue",
            "1",
            "--out",
            explicit_out.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!explicit.status.success());
    let stderr = String::from_utf8_lossy(&explicit.stderr);
    assert!(stderr.contains("GFF has no gene_biotype"), "{stderr}");
    assert!(stderr.contains("provide a GFF with gene_biotype"), "{stderr}");

    let pack_out = dir.join("pack-default.tsv");
    let pack_report = dir.join("pack-default-report.json");
    let pack = Command::new(bin)
        .args([
            "scan",
            "--genome",
            genome_path.to_str().unwrap(),
            "--gff",
            gff_path.to_str().unwrap(),
            "--motifs",
            motif_path.to_str().unwrap(),
            "--gene-biotype",
            "protein_coding",
            "--gene-biotype-from-pack",
            "--promoter-len",
            "100",
            "--pvalue",
            "1",
            "--out",
            pack_out.to_str().unwrap(),
            "--report-json",
            pack_report.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        pack.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&pack.stderr)
    );
    let text = std::fs::read_to_string(&pack_out).unwrap();
    assert!(text.contains("geneA"), "{text}");
    assert!(text.contains("geneB"), "{text}");
    let report: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&pack_report).unwrap()).unwrap();
    assert!(
        report["diagnostics"]
            .as_array()
            .unwrap()
            .iter()
            .any(|d| d["code"] == "GENE_BIOTYPE_FILTER_SKIPPED"),
        "{report}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn cli_trace_outputs() {
    let dir = std::env::temp_dir().join(format!(
        "regutrace_trace_{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let genes = dir.join("tf_genes.tsv");
    write(
        &genes,
        "gene_id\ttf_symbol\tfamily\tidentification_status\tsource\n\
TF1\tTF1\tF\tidentified\ttest\n\
TF2\tTF2\tF\tidentified\ttest\n",
    );
    let map = dir.join("tf_map.tsv");
    write(
        &map,
        "tf_id\ttf_symbol\tfamily\tmotif_id\treference_tf_id\tmapping_method\tdomain_mode\textraction_family\textraction_models\tdomain_count\n\
TF1\tTF1\tF\tMA1\tR1\treciprocal_dbd_similarity\tfull_architecture\tF\tMODEL\t1\n\
TF2\tTF2\tF\tMA1\tR1\treciprocal_dbd_tied_candidate\tfull_architecture\tF\tMODEL\t1\n",
    );
    let reference = dir.join("reference-map.tsv");
    write(
        &reference,
        "tf_id\ttf_symbol\ttf_family\tmotif_id\n\
R1\tREF1\tF\tMA1\n",
    );
    let motifs = dir.join("motifs.txt");
    write(
        &motifs,
        ">MA1 MOTIF ONE\nA [ 10 0 ]\nC [ 0 10 ]\nG [ 0 0 ]\nT [ 0 0 ]\n",
    );
    let evidence = dir.join("evidence.tsv");
    write(
        &evidence,
        "tf_id\ttf_symbol\ttarget_id\ttarget_symbol\tevidence_type\tevidence_level\tmethod\tsource\tpmid\tcondition\tpromoter_start\tpromoter_end\tmotif_id\thit_start\thit_end\tstrand\tscore\tp_value\tnotes\n\
TF1\tTF1\tG1\tGENE1\tprediction_motif\tD\tjaspar_pwm\tpromoter_scan\t\t\t100\t200\tMA1\t10\t20\t+\t10\t1e-5\t\n\
TF2\tTF2\tG1\tGENE1\tprediction_motif\tD\tjaspar_pwm\tpromoter_scan\t\t\t100\t200\tMA1\t10\t20\t+\t10\t1e-5\t\n\
TF1\tTF1\tG2\tGENE2\tprediction_motif\tD\tjaspar_pwm\tpromoter_scan\t\t\t300\t400\tMA1\t30\t40\t+\t8\t1e-4\t\n",
    );
    let out = dir.join("trace");
    let bin = env!("CARGO_BIN_EXE_regutrace");
    let result = Command::new(bin)
        .args([
            "trace",
            "--tf-genes",
            genes.to_str().unwrap(),
            "--tf-map",
            map.to_str().unwrap(),
            "--evidence",
            evidence.to_str().unwrap(),
            "--motifs",
            motifs.to_str().unwrap(),
            "--reference-map",
            reference.to_str().unwrap(),
            "--pack-version",
            "test-pack",
            "--out-dir",
            out.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&result.stderr)
    );

    let trace = std::fs::read_to_string(out.join("tf_motif_trace.tsv")).unwrap();
    assert!(trace.contains("MOTIF ONE"), "{trace}");
    assert!(trace.contains("REF1"), "{trace}");
    assert!(trace.contains("test-pack"), "{trace}");
    let tf2 = trace.lines().find(|line| line.starts_with("TF2\t")).unwrap();
    assert!(tf2.contains("reciprocal_dbd_tied_candidate"), "{tf2}");
    assert!(tf2.contains("\ttrue\t"), "{tf2}");
    assert!(tf2.contains("\t2\t"), "{tf2}");

    let ambiguity = std::fs::read_to_string(out.join("motif_ambiguity.tsv")).unwrap();
    assert!(ambiguity.contains("MA1"), "{ambiguity}");
    assert!(ambiguity.contains("TF1;TF2"), "{ambiguity}");
    assert!(ambiguity.contains("REF1"), "{ambiguity}");

    let ties = std::fs::read_to_string(out.join("target_candidate_ties.tsv")).unwrap();
    assert!(ties.contains("G1"), "{ties}");
    assert!(ties.contains("TF1;TF2"), "{ties}");
    assert!(ties.contains("GENE1"), "{ties}");

    let _ = std::fs::remove_dir_all(&dir);
}
#[test]
fn phase0_mapping_diagnostics_are_auditable() {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/phase0");
    let dir = std::env::temp_dir().join(format!(
        "regutrace_phase0_{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let map = dir.join("tf_map.tsv");
    let diagnostics = dir.join("mapping-diagnostics.tsv");
    let families = dir.join("mapping-diagnostics.family.tsv");
    let bin = env!("CARGO_BIN_EXE_regutrace");
    let result = Command::new(bin)
        .args([
            "map-tf-motifs",
            "--tf-genes",
            fixture.join("tf_genes.tsv").to_str().unwrap(),
            "--reference-map",
            fixture.join("reference-map.tsv").to_str().unwrap(),
            "--forward",
            fixture.join("forward.tsv").to_str().unwrap(),
            "--reverse",
            fixture.join("reverse.tsv").to_str().unwrap(),
            "--out",
            map.to_str().unwrap(),
            "--diagnostics-out",
            diagnostics.to_str().unwrap(),
            "--diagnostics-family-out",
            families.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&result.stderr)
    );

    let diagnostics = std::fs::read_to_string(diagnostics).unwrap();
    let row = |id: &str| {
        diagnostics
            .lines()
            .find(|line| line.starts_with(&format!("{id}\t")))
            .unwrap_or_else(|| panic!("missing diagnostic row {id}"))
            .split('\t')
            .collect::<Vec<_>>()
    };
    assert_eq!(row("G_MAPPED")[5], "MAPPED");
    assert_eq!(row("G_TIED")[5], "MAPPED_TIED");
    assert_eq!(row("G_NODONOR")[5], "NO_REFERENCE_DONOR");
    assert_eq!(row("G_NOBLAST")[5], "NO_BLAST_HIT");
    assert_eq!(row("G_BELOW")[5], "BELOW_IDENTITY_OR_COVERAGE");
    assert_eq!(row("G_NONRECIP")[5], "NOT_RECIPROCAL_BEST");
    assert_eq!(row("G_FAMILY")[5], "NO_REFERENCE_DONOR");
    assert!(row("G_FAMILY")[6].contains("FAMILY_MISMATCH"));
    assert_eq!(row("G_AMB")[5], "AMBIGUOUS_EXCLUDED");
    assert_eq!(row("G_REJ")[5], "REJECTED_EXCLUDED");

    let family_text = std::fs::read_to_string(families).unwrap();
    let family = |name: &str| {
        family_text
            .lines()
            .find(|line| line.starts_with(&format!("{name}\t")))
            .unwrap_or_else(|| panic!("missing family row {name}"))
            .split('\t')
            .collect::<Vec<_>>()
    };
    let a = family("A");
    assert_eq!(a[1], "7");
    assert_eq!(a[2], "5");
    assert_eq!(a[5], "2");
    assert_eq!(a[7], "3");
    let reason_total: usize = a[10..].iter().map(|v| v.parse::<usize>().unwrap()).sum();
    assert_eq!(reason_total, 3, "{a:?}");
    assert_eq!(family("B")[10], "1");
    assert_eq!(family("C")[10], "1");
    assert!(!map.as_os_str().is_empty());
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn registry_inspect_accepts_legacy_pack_and_writes_summary() {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/phase0");
    let dir = std::env::temp_dir().join(format!(
        "regutrace_registry_{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::copy(fixture.join("reference-map.tsv"), dir.join("reference-map.tsv")).unwrap();
    write(
        &dir.join("manifest.json"),
        "{\"pack_version\":\"test-pack-1\",\"species\":\"synthetic\",\
         \"files\":{\"reference_map\":\"reference-map.tsv\"}}",
    );
    write(
        &dir.join("motifs.txt"),
        ">MA1 a\nA [ 1 0 ]\nC [ 0 1 ]\nG [ 0 0 ]\nT [ 0 0 ]\n\
>MA2 b\nA [ 1 0 ]\nC [ 0 1 ]\nG [ 0 0 ]\nT [ 0 0 ]\n\
>MA3 c\nA [ 0 1 ]\nC [ 0 1 ]\nG [ 1 0 ]\nT [ 1 0 ]\n\
>MA4 d\nA [ 1 1 ]\nC [ 0 0 ]\nG [ 0 0 ]\nT [ 0 0 ]\n",
    );
    let registry = dir.join("registry.tsv");
    let summary = dir.join("summary.json");
    let bin = env!("CARGO_BIN_EXE_regutrace");
    let result = Command::new(bin)
        .args([
            "registry",
            "inspect",
            "--reference-pack",
            dir.to_str().unwrap(),
            "--out",
            registry.to_str().unwrap(),
            "--summary-json",
            summary.to_str().unwrap(),
            "--motif-library",
            dir.join("motifs.txt").to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    let text = std::fs::read_to_string(&registry).unwrap();
    let rows: Vec<_> = text.lines().collect();
    assert_eq!(rows.len(), 4, "{text}");
    assert!(rows[0].starts_with("motif_id\tdonor_id\tfamily\tsource"));
    assert!(rows[1].contains("\tD1\tA\tlegacy-pack\ttest-pack-1\t"), "{}", rows[1]);
    let summary: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&summary).unwrap()).unwrap();
    assert_eq!(summary["records"], 3);
    assert_eq!(summary["unique_motifs"], 3);
    assert_eq!(summary["unique_donors"], 3);
    assert_eq!(summary["missing_dbd"], 3);
    assert_eq!(summary["qc_status"]["legacy_unverified"], 3);
    assert_eq!(summary["donor_species"]["synthetic"], 3);
    assert_eq!(summary["families"][0]["family"], "A");
    assert_eq!(summary["library"]["motifs"], 4);
    assert_eq!(summary["library"]["library_unmapped"], 1);
    assert_eq!(summary["library"]["duplicate_group_count"], 1);
    assert_eq!(summary["library"]["duplicate_motif_count"], 2);
    let unknown = Command::new(bin)
        .args(["registry", "bogus", "--out", registry.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!unknown.status.success());
    assert!(String::from_utf8_lossy(&unknown.stderr).contains("unknown registry command"));
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn phase2_dbd_extract_and_mapping_are_auditable() {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/phase2");
    let dir = std::env::temp_dir().join(format!(
        "regutrace_phase2_{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let bin = env!("CARGO_BIN_EXE_regutrace");
    let local_faa = dir.join("local-dbd.faa");
    let local_tsv = dir.join("local-dbd.tsv");
    let extracted = Command::new(bin)
        .args([
            "dbd-extract",
            "--tf-genes",
            fixture.join("tf_genes.tsv").to_str().unwrap(),
            "--proteins",
            fixture.join("proteins.faa").to_str().unwrap(),
            "--dbd-models",
            fixture.join("dbd-models.json").to_str().unwrap(),
            "--out-faa",
            local_faa.to_str().unwrap(),
            "--out-tsv",
            local_tsv.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        extracted.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&extracted.stderr)
    );
    let fasta = std::fs::read_to_string(&local_faa).unwrap();
    assert_eq!(fasta.matches('>').count(), 4, "{fasta}");
    assert!(fasta.contains(">G1\nACGTACGTAC\n"), "{fasta}");
    let table = std::fs::read_to_string(&local_tsv).unwrap();
    assert_eq!(table.lines().count(), 5);
    let g1 = table
        .lines()
        .find(|line| line.starts_with("G1\t"))
        .unwrap()
        .split('\t')
        .collect::<Vec<_>>();
    assert_eq!(g1[4], "3");
    assert_eq!(g1[5], "12");
    assert_eq!(g1[6], "1");
    assert_eq!(g1[7], "full_architecture");

    let map = dir.join("dbd-map.tsv");
    let diagnostics = dir.join("dbd-diagnostics.tsv");
    let mapped = Command::new(bin)
        .args([
            "map-dbd-motifs",
            "--tf-genes",
            fixture.join("tf_genes.tsv").to_str().unwrap(),
            "--proteins",
            fixture.join("proteins.faa").to_str().unwrap(),
            "--dbd-models",
            fixture.join("dbd-models.json").to_str().unwrap(),
            "--reference-map",
            fixture.join("reference-map.tsv").to_str().unwrap(),
            "--donor-dbd",
            fixture.join("donor-dbd.tsv").to_str().unwrap(),
            "--donor-dbd-faa",
            fixture.join("donor-dbd.faa").to_str().unwrap(),
            "--threshold-config",
            fixture.join("thresholds.json").to_str().unwrap(),
            "--allow-uncalibrated",
            "true",
            "--dbd-alignments",
            fixture.join("dbd-alignments.tsv").to_str().unwrap(),
            "--out",
            map.to_str().unwrap(),
            "--diagnostics-out",
            diagnostics.to_str().unwrap(),
            "--reference-pack-version",
            "test-pack",
        ])
        .output()
        .unwrap();
    assert!(
        mapped.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&mapped.stderr)
    );
    let map_text = std::fs::read_to_string(&map).unwrap();
    let rows = map_text
        .lines()
        .skip(1)
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .collect::<Vec<_>>();
    assert_eq!(rows.len(), 3, "{map_text}");
    assert_eq!(rows.iter().filter(|row| row[0] == "G1").count(), 2);
    assert_eq!(rows.iter().filter(|row| row[0] == "G3").count(), 1);
    let g1_rows = rows.iter().filter(|row| row[0] == "G1").collect::<Vec<_>>();
    assert!(g1_rows.iter().all(|row| row[11] == "dbd_exploratory"));
    assert!(g1_rows.iter().all(|row| row[32] == "multiple_donors_and_motifs"));
    assert!(g1_rows.iter().all(|row| row[33] == "DBD_TRANSFER_AMBIGUOUS"));
    assert!(g1_rows.iter().all(|row| row[34] == "test-pack"));
    let g3 = rows.iter().find(|row| row[0] == "G3").unwrap();
    assert_eq!(g3[11], "dbd_exploratory");
    assert_eq!(g3[32], "unique");
    assert_eq!(g3[33], "DBD_TRANSFER");

    let diagnostic_text = std::fs::read_to_string(&diagnostics).unwrap();
    let reason = |id: &str| {
        diagnostic_text
            .lines()
            .find(|line| line.starts_with(&format!("{id}\t")))
            .unwrap_or_else(|| panic!("missing diagnostic for {id}"))
            .split('\t')
            .collect::<Vec<_>>()[5]
            .to_owned()
    };
    assert_eq!(reason("G1"), "DBD_EXPLORATORY");
    assert_eq!(reason("G2"), "BELOW_DBD_THRESHOLD");
    assert_eq!(reason("G3"), "DBD_EXPLORATORY");
    assert_eq!(reason("G4"), "NO_DBD_SEQUENCE");
    assert_eq!(reason("G5"), "NO_BLAST_HIT");
    let g1_diag = diagnostic_text
        .lines()
        .find(|line| line.starts_with("G1\t"))
        .unwrap()
        .split('\t')
        .collect::<Vec<_>>();
    assert_eq!(g1_diag[12], "2");
    assert_eq!(g1_diag[13], "D1");
    assert_eq!(g1_diag[14], "90.0000");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn phase2_5_calibrated_mapping_uses_reciprocity_and_specificity_groups() {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/phase2");
    let dir = std::env::temp_dir().join(format!(
        "regutrace_phase2_5_{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let bin = env!("CARGO_BIN_EXE_regutrace");
    let map = dir.join("calibrated-map.tsv");
    let diagnostics = dir.join("calibrated-diagnostics.tsv");
    let output = Command::new(bin)
        .args([
            "map-dbd-motifs",
            "--tf-genes",
            fixture.join("tf_genes.tsv").to_str().unwrap(),
            "--proteins",
            fixture.join("proteins.faa").to_str().unwrap(),
            "--dbd-models",
            fixture.join("dbd-models.json").to_str().unwrap(),
            "--reference-map",
            fixture.join("reference-map.tsv").to_str().unwrap(),
            "--donor-dbd",
            fixture.join("donor-dbd.tsv").to_str().unwrap(),
            "--donor-dbd-faa",
            fixture.join("donor-dbd.faa").to_str().unwrap(),
            "--threshold-config",
            fixture.join("thresholds-calibrated.json").to_str().unwrap(),
            "--specificity-groups",
            fixture.join("specificity-groups.tsv").to_str().unwrap(),
            "--require-reciprocal-dbd",
            "true",
            "--dbd-alignments",
            fixture.join("dbd-alignments.tsv").to_str().unwrap(),
            "--reverse-dbd-alignments",
            fixture.join("reverse-dbd-alignments.tsv").to_str().unwrap(),
            "--out",
            map.to_str().unwrap(),
            "--diagnostics-out",
            diagnostics.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let map_text = std::fs::read_to_string(&map).unwrap();
    let rows = map_text
        .lines()
        .skip(1)
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .collect::<Vec<_>>();
    assert_eq!(rows.len(), 3);
    assert!(rows.iter().all(|row| row[10] == "calibrated"));
    assert!(rows.iter().all(|row| row[11] == "dbd_calibrated"));
    assert!(rows.iter().all(|row| row[12] == "yes"));
    let g1 = rows.iter().filter(|row| row[0] == "G1").collect::<Vec<_>>();
    assert_eq!(g1.len(), 2);
    assert!(g1.iter().all(|row| row[13] == "A1"));
    assert!(g1.iter().all(|row| row[14] == "1"));
    let g3 = rows.iter().find(|row| row[0] == "G3").unwrap();
    assert_eq!(g3[13], "B1");
    assert_eq!(g3[33], "DBD_CALIBRATED");
    let diagnostic_text = std::fs::read_to_string(&diagnostics).unwrap();
    let reason = |id: &str| {
        diagnostic_text
            .lines()
            .find(|line| line.starts_with(&format!("{id}\t")))
            .unwrap()
            .split('\t')
            .collect::<Vec<_>>()[5]
            .to_owned()
    };
    assert_eq!(reason("G1"), "DBD_CALIBRATED");
    assert_eq!(reason("G3"), "DBD_CALIBRATED");

    let missing_groups = Command::new(bin)
        .args([
            "map-dbd-motifs",
            "--tf-genes",
            fixture.join("tf_genes.tsv").to_str().unwrap(),
            "--proteins",
            fixture.join("proteins.faa").to_str().unwrap(),
            "--dbd-models",
            fixture.join("dbd-models.json").to_str().unwrap(),
            "--reference-map",
            fixture.join("reference-map.tsv").to_str().unwrap(),
            "--donor-dbd",
            fixture.join("donor-dbd.tsv").to_str().unwrap(),
            "--donor-dbd-faa",
            fixture.join("donor-dbd.faa").to_str().unwrap(),
            "--threshold-config",
            fixture.join("thresholds-calibrated.json").to_str().unwrap(),
            "--dbd-alignments",
            fixture.join("dbd-alignments.tsv").to_str().unwrap(),
            "--out",
            dir.join("should-not-exist.tsv").to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!missing_groups.status.success());
    assert!(
        String::from_utf8_lossy(&missing_groups.stderr).contains("requires --specificity-groups")
    );
    let missing_reciprocal = Command::new(bin)
        .args([
            "map-dbd-motifs",
            "--tf-genes",
            fixture.join("tf_genes.tsv").to_str().unwrap(),
            "--proteins",
            fixture.join("proteins.faa").to_str().unwrap(),
            "--dbd-models",
            fixture.join("dbd-models.json").to_str().unwrap(),
            "--reference-map",
            fixture.join("reference-map.tsv").to_str().unwrap(),
            "--donor-dbd",
            fixture.join("donor-dbd.tsv").to_str().unwrap(),
            "--donor-dbd-faa",
            fixture.join("donor-dbd.faa").to_str().unwrap(),
            "--threshold-config",
            fixture.join("thresholds-calibrated.json").to_str().unwrap(),
            "--specificity-groups",
            fixture.join("specificity-groups.tsv").to_str().unwrap(),
            "--dbd-alignments",
            fixture.join("dbd-alignments.tsv").to_str().unwrap(),
            "--out",
            dir.join("should-not-exist-2.tsv").to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!missing_reciprocal.status.success());
    assert!(
        String::from_utf8_lossy(&missing_reciprocal.stderr)
            .contains("requires --require-reciprocal-dbd")
    );
    let _ = std::fs::remove_dir_all(&dir);
}
