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
            "--pvalue",
            "1e-3",
            "--annotation",
            anno_path.to_str().unwrap(),
            "--user-list",
            user_path.to_str().unwrap(),
            "--out",
            out_path.to_str().unwrap(),
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
