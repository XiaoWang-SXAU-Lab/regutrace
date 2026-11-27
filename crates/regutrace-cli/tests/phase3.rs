//! Phase 3 scan architecture tests: exact PWM reuse, a site layer separate
//! from TF evidence rows, persistent cache invalidation, and empirical
//! background provenance.

use std::path::PathBuf;
use std::process::Command;

fn write(path: &PathBuf, text: &str) {
    std::fs::write(path, text).unwrap();
}

fn pwm_library(ids: &[&str], motif: &str) -> String {
    let mut text = String::new();
    for id in ids {
        text.push_str(&format!(">{id} TEST\n"));
        for base in ["A", "C", "G", "T"] {
            let values: Vec<String> = motif
                .chars()
                .map(|letter| {
                    if letter.to_string() == base {
                        "10".to_string()
                    } else {
                        "0".to_string()
                    }
                })
                .collect();
            text.push_str(&format!("{base} [ {} ]\n", values.join(" ")));
        }
    }
    text
}

#[test]
fn phase3_exact_pwm_reuse_sites_cache_and_empirical_background() {
    let dir = std::env::temp_dir().join(format!("regutrace_phase3_cli_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let motif = "ACGTAC";
    let mut genome = vec![b'A'; 2000];
    genome[150..156].copy_from_slice(motif.as_bytes());
    let genome_path = dir.join("genome.fa");
    write(
        &genome_path,
        &format!(">Chr1\n{}\n", String::from_utf8(genome).unwrap()),
    );
    let gff_path = dir.join("annotation.gff3");
    write(
        &gff_path,
        "##gff-version 3\n\
Chr1\ttest\tgene\t201\t400\t.\t+\t.\tID=gene:G1;Name=GENE1\n\
Chr1\ttest\tmRNA\t201\t400\t.\t+\t.\tID=transcript:G1.t1;Parent=gene:G1\n",
    );
    let motif_path = dir.join("motifs.txt");
    write(&motif_path, &pwm_library(&["MA0001.1", "MA0001.2"], motif));
    let map_path = dir.join("tf_map.tsv");
    write(
        &map_path,
        "tf_id\ttf_symbol\tfamily\tmotif_id\n\
TF1\tTF1\tTEST\tMA0001.1\n\
TF2\tTF2\tTEST\tMA0001.2\n",
    );
    let cache_dir = dir.join("cache");
    let bin = env!("CARGO_BIN_EXE_regutrace");

    let run = |name: &str, extra: &[&str]| {
        let out = dir.join(format!("{name}.tsv"));
        let sites = dir.join(format!("{name}-sites.tsv"));
        let testability = dir.join(format!("{name}-testability.tsv"));
        let summary = dir.join(format!("{name}-promoters.tsv"));
        let manifest = dir.join(format!("{name}-manifest.json"));
        let mut args = vec![
            "scan".to_string(),
            "--genome".to_string(),
            genome_path.to_string_lossy().into_owned(),
            "--gff".to_string(),
            gff_path.to_string_lossy().into_owned(),
            "--motifs".to_string(),
            motif_path.to_string_lossy().into_owned(),
            "--tf-map".to_string(),
            map_path.to_string_lossy().into_owned(),
            "--promoter-len".to_string(),
            "100".to_string(),
            "--downstream".to_string(),
            "100".to_string(),
            "--pvalue".to_string(),
            "1e-3".to_string(),
            "--reuse-exact-pwm".to_string(),
            "--scan-cache".to_string(),
            cache_dir.to_string_lossy().into_owned(),
            "--sites-out".to_string(),
            sites.to_string_lossy().into_owned(),
            "--motif-testability-out".to_string(),
            testability.to_string_lossy().into_owned(),
            "--promoter-summary-out".to_string(),
            summary.to_string_lossy().into_owned(),
            "--scan-manifest".to_string(),
            manifest.to_string_lossy().into_owned(),
            "--out".to_string(),
            out.to_string_lossy().into_owned(),
        ];
        args.extend(extra.iter().map(|value| value.to_string()));
        let output = Command::new(bin).args(&args).output().unwrap();
        assert!(
            output.status.success(),
            "stderr: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        (
            out,
            sites,
            testability,
            summary,
            manifest,
            String::from_utf8_lossy(&output.stderr).into_owned(),
        )
    };

    // The block below re-baselines the legacy uniform scan: its site rows,
    // evidence counts and cache-key expectations were all measured under that
    // null, so it is pinned explicitly now that the default is the fitted
    // promoter composition.
    let (out, sites, testability, summary, manifest, stderr) =
        run("first", &["--background-mode", "uniform"]);
    let evidence = std::fs::read_to_string(&out).unwrap();
    assert_eq!(evidence.lines().count().saturating_sub(1), 2, "{stderr}");
    let site_text = std::fs::read_to_string(&sites).unwrap();
    assert_eq!(site_text.lines().count().saturating_sub(1), 1, "{site_text}");
    let site_row: Vec<&str> = site_text.lines().nth(1).unwrap().split('\t').collect();
    assert_eq!(site_row[1], "MA0001.1;MA0001.2");
    assert_eq!(site_row[10], "151");
    assert_eq!(site_row[11], "156");
    assert_eq!(site_row[12], "+");
    let summary_text = std::fs::read_to_string(&summary).unwrap();
    let summary_row: Vec<&str> = summary_text.lines().nth(1).unwrap().split('\t').collect();
    assert_eq!(summary_row[7], "1");
    assert_eq!(summary_row[9], "1");
    let testability_text = std::fs::read_to_string(&testability).unwrap();
    assert_eq!(testability_text.lines().count().saturating_sub(1), 2);
    let manifest_json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&manifest).unwrap()).unwrap();
    assert_eq!(manifest_json["motifs"]["pwm_groups"], 1);
    assert_eq!(manifest_json["motifs"]["unique_pwms"], 1);
    assert_eq!(manifest_json["motifs"]["duplicate_pwm_groups"], 1);
    assert_eq!(manifest_json["sites"]["unique_pwm_sites"], 1);
    assert_eq!(manifest_json["sites"]["evidence_rows"], 2);
    assert_eq!(manifest_json["scan_jobs"]["cache_misses"], 1);
    assert_eq!(manifest_json["scan_jobs"]["cache_hits"], 0);

    // The same cache and inputs must be a pure cache hit.
    let (_, _, _, _, second_manifest, _) =
        run("second", &["--background-mode", "uniform"]);
    let second: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&second_manifest).unwrap()).unwrap();
    assert_eq!(second["scan_jobs"]["cache_hits"], 1);
    assert_eq!(second["scan_jobs"]["cache_misses"], 0);

    // Changing the requested p-value must change the cache key.
    let (_, _, _, _, changed_manifest, _) = run(
        "changed",
        &["--background-mode", "uniform", "--pvalue", "1e-2"],
    );
    let changed: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&changed_manifest).unwrap()).unwrap();
    assert_eq!(changed["scan_jobs"]["cache_hits"], 0);
    assert_eq!(changed["scan_jobs"]["cache_misses"], 1);

    // The default null is the fitted promoter composition. Naming it explicitly
    // differs from the pinned uniform runs above, and the manifest must record
    // which null produced the p-values.
    let (_, _, _, _, empirical_manifest, _) = run(
        "empirical",
        &["--background-mode", "empirical", "--pvalue", "1e-3"],
    );
    let empirical: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&empirical_manifest).unwrap()).unwrap();
    assert_eq!(empirical["scan_jobs"]["cache_hits"], 0);
    assert_eq!(empirical["scan_jobs"]["cache_misses"], 1);
    assert!(empirical["background"]["mode"]
        .as_str()
        .unwrap()
        .starts_with("empirical"));
    assert_eq!(empirical["null_model"]["calibrated"], true);
    assert_eq!(
        empirical["background"]["calibrated_to_promoter_composition"],
        true
    );
    let bg: Vec<f64> = empirical["background"]["values_acgt"]
        .as_array()
        .unwrap()
        .iter()
        .map(|value| value.as_f64().unwrap())
        .collect();
    assert!((bg.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    assert_eq!(bg[0], bg[3]);
    assert_eq!(bg[1], bg[2]);

    // Omitting the mode must land on the same cache entry as naming empirical:
    // the default is that model, not a third one.
    let (_, _, _, _, default_manifest, _) = run("default-null", &["--pvalue", "1e-3"]);
    let default_null: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&default_manifest).unwrap()).unwrap();
    assert!(default_null["background"]["mode"]
        .as_str()
        .unwrap()
        .starts_with("empirical"));
    assert_eq!(default_null["scan_jobs"]["cache_hits"], 1);
    assert_eq!(default_null["scan_jobs"]["cache_misses"], 0);

    // A 2-bp motif cannot reach p <= 1e-8 and must be reported, not silently
    // treated as a motif with zero biological signal.
    let short_motif_path = dir.join("short-motifs.txt");
    write(
        &short_motif_path,
        ">MA9999.1 SHORT\nA [ 10 0 ]\nC [ 0 10 ]\nG [ 0 0 ]\nT [ 0 0 ]\n",
    );
    let short_map_path = dir.join("short-tf_map.tsv");
    write(
        &short_map_path,
        "tf_id\ttf_symbol\tfamily\tmotif_id\nSHORT_TF\tSHORT_TF\tTEST\tMA9999.1\n",
    );
    let unattainable_out = dir.join("unattainable.tsv");
    let unattainable_testability = dir.join("unattainable-testability.tsv");
    let output = Command::new(bin)
        .args([
            "scan",
            "--genome",
            genome_path.to_str().unwrap(),
            "--gff",
            gff_path.to_str().unwrap(),
            "--motifs",
            short_motif_path.to_str().unwrap(),
            "--tf-map",
            short_map_path.to_str().unwrap(),
            "--promoter-len",
            "100",
            "--downstream",
            "100",
            "--pvalue",
            "1e-8",
            "--reuse-exact-pwm",
            "--motif-testability-out",
            unattainable_testability.to_str().unwrap(),
            "--out",
            unattainable_out.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let table = std::fs::read_to_string(&unattainable_testability).unwrap();
    let data_row: Vec<&str> = table.lines().nth(1).unwrap().split('\t').collect();
    assert_eq!(data_row[8], "threshold_unattainable");
    assert_eq!(data_row[7], "false");
    let evidence = std::fs::read_to_string(&unattainable_out).unwrap();
    assert_eq!(evidence.lines().count(), 1);
}

#[test]
fn enrichment_interface_uses_scanned_universe_and_reports_bh() {
    let dir = std::env::temp_dir().join(format!(
        "regutrace_phase3_enrich_{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let sites = dir.join("sites.tsv");
    let mut site_text = "pwm_hash\tmotif_ids\tgene_id\n".to_string();
    for gene in 1..=5 {
        site_text.push_str(&format!("h1\tM1\tG{gene}\n"));
    }
    site_text.push_str("h2\tM2\tG1\n");
    write(&sites, &site_text);
    let summary = dir.join("promoter-summary.tsv");
    let mut summary_text = "gene_id\n".to_string();
    for gene in 1..=10 {
        summary_text.push_str(&format!("G{gene}\n"));
    }
    write(&summary, &summary_text);
    let gmt = dir.join("sets.gmt");
    write(&gmt, "SET_A\tsynthetic\tG1\tG2\tG3\tG4\tG5\n");
    let out = dir.join("enrichment.tsv");
    let bin = env!("CARGO_BIN_EXE_regutrace");
    let output = Command::new(bin)
        .args([
            "enrich",
            "--sites",
            sites.to_str().unwrap(),
            "--promoter-summary",
            summary.to_str().unwrap(),
            "--gene-sets",
            gmt.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(&out).unwrap();
    let rows: Vec<Vec<&str>> = text
        .lines()
        .skip(1)
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), 2, "{text}");
    let m1 = rows.iter().find(|row| row[4] == "M1").unwrap();
    let m2 = rows.iter().find(|row| row[4] == "M2").unwrap();
    assert_eq!(m1[6], "5");
    assert_eq!(m2[6], "1");
    assert!(m1[11].parse::<f64>().unwrap() <= 0.05, "{text}");
    assert!(m2[11].parse::<f64>().unwrap() > 0.05, "{text}");
    let _ = std::fs::remove_dir_all(&dir);
}
