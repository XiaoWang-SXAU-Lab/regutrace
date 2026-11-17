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
