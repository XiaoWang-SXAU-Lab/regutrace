//! SQLite evidence store. One row per evidence record; duplicate imports are
//! ignored via a stable record key.

use std::path::Path;

use anyhow::{Context, Result};
use rusqlite::{params, Connection};

use crate::io::EdgeRecord;

const SCHEMA: &str = r#"
CREATE TABLE IF NOT EXISTS evidence (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    record_key TEXT NOT NULL UNIQUE,
    tf_id TEXT NOT NULL,
    tf_symbol TEXT NOT NULL DEFAULT '',
    target_id TEXT NOT NULL,
    target_symbol TEXT NOT NULL DEFAULT '',
    evidence_type TEXT NOT NULL,
    evidence_level TEXT NOT NULL,
    method TEXT NOT NULL DEFAULT '',
    source TEXT NOT NULL DEFAULT '',
    pmid TEXT NOT NULL DEFAULT '',
    condition TEXT NOT NULL DEFAULT '',
    promoter_start INTEGER,
    promoter_end INTEGER,
    motif_id TEXT NOT NULL DEFAULT '',
    hit_start INTEGER,
    hit_end INTEGER,
    strand TEXT NOT NULL DEFAULT '',
    score REAL,
    p_value REAL,
    notes TEXT NOT NULL DEFAULT '',
    imported_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
);
CREATE INDEX IF NOT EXISTS idx_evidence_tf ON evidence(tf_id);
CREATE INDEX IF NOT EXISTS idx_evidence_target ON evidence(target_id);
CREATE INDEX IF NOT EXISTS idx_evidence_level ON evidence(evidence_level);
"#;

pub struct Store {
    conn: Connection,
}

fn opt_u64(v: Option<u64>) -> String {
    v.map(|x| x.to_string()).unwrap_or_default()
}
fn record_key(e: &EdgeRecord) -> String {
    format!(
        "{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
        e.tf_id,
        e.target_id,
        e.evidence_type,
        e.method,
        e.source,
        e.motif_id,
        opt_u64(e.promoter_start),
        opt_u64(e.hit_start),
        e.strand,
        e.condition
    )
}

impl Store {
    pub fn open(path: &Path) -> Result<Self> {
        let conn = Connection::open(path)
            .with_context(|| format!("open evidence database {}", path.display()))?;
        conn.execute_batch(SCHEMA)?;
        Ok(Self { conn })
    }

    /// Insert evidence rows; returns the number of newly inserted rows.
    pub fn insert_edges(&mut self, edges: &[EdgeRecord]) -> Result<usize> {
        let tx = self.conn.transaction()?;
        let mut inserted = 0usize;
        {
            let mut stmt = tx.prepare(
                "INSERT OR IGNORE INTO evidence (
                    record_key, tf_id, tf_symbol, target_id, target_symbol,
                    evidence_type, evidence_level, method, source, pmid, condition,
                    promoter_start, promoter_end, motif_id, hit_start, hit_end, strand,
                    score, p_value, notes
                ) VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11,?12,?13,?14,?15,?16,?17,?18,?19,?20)",
            )?;
            for e in edges {
                inserted += stmt.execute(params![
                    record_key(e),
                    e.tf_id,
                    e.tf_symbol,
                    e.target_id,
                    e.target_symbol,
                    e.evidence_type,
                    e.evidence_level,
                    e.method,
                    e.source,
                    e.pmid,
                    e.condition,
                    e.promoter_start.map(|v| v as i64),
                    e.promoter_end.map(|v| v as i64),
                    e.motif_id,
                    e.hit_start.map(|v| v as i64),
                    e.hit_end.map(|v| v as i64),
                    e.strand,
                    e.score,
                    e.p_value,
                    e.notes,
                ])?;
            }
        }
        tx.commit()?;
        Ok(inserted)
    }

    pub fn export_edges(&self) -> Result<Vec<EdgeRecord>> {
        let mut stmt = self.conn.prepare(
            "SELECT tf_id, tf_symbol, target_id, target_symbol, evidence_type,
                    evidence_level, method, source, pmid, condition, promoter_start,
                    promoter_end, motif_id, hit_start, hit_end, strand, score, p_value, notes
             FROM evidence ORDER BY id",
        )?;
        let rows = stmt.query_map([], |row| {
            let opt_i64 = |i: usize| row.get::<_, Option<i64>>(i);
            Ok(EdgeRecord {
                tf_id: row.get(0)?,
                tf_symbol: row.get(1)?,
                target_id: row.get(2)?,
                target_symbol: row.get(3)?,
                evidence_type: row.get(4)?,
                evidence_level: row.get(5)?,
                method: row.get(6)?,
                source: row.get(7)?,
                pmid: row.get(8)?,
                condition: row.get(9)?,
                promoter_start: opt_i64(10)?.map(|v| v as u64),
                promoter_end: opt_i64(11)?.map(|v| v as u64),
                motif_id: row.get(12)?,
                hit_start: opt_i64(13)?.map(|v| v as u64),
                hit_end: opt_i64(14)?.map(|v| v as u64),
                strand: row.get(15)?,
                score: row.get(16)?,
                p_value: row.get(17)?,
                notes: row.get(18)?,
            })
        })?;
        let mut out = Vec::new();
        for row in rows {
            out.push(row?);
        }
        Ok(out)
    }

    /// (evidence_level, evidence_type, count)
    pub fn stats(&self) -> Result<Vec<(String, String, i64)>> {
        let mut stmt = self.conn.prepare(
            "SELECT evidence_level, evidence_type, COUNT(*)
             FROM evidence GROUP BY evidence_level, evidence_type
             ORDER BY evidence_level, evidence_type",
        )?;
        let rows = stmt.query_map([], |row| {
            Ok((row.get(0)?, row.get(1)?, row.get(2)?))
        })?;
        let mut out = Vec::new();
        for row in rows {
            out.push(row?);
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn edge() -> EdgeRecord {
        EdgeRecord {
            tf_id: "TF1".into(),
            tf_symbol: "TF1".into(),
            target_id: "G1".into(),
            target_symbol: "GENE1".into(),
            evidence_type: "prediction_motif".into(),
            evidence_level: "D".into(),
            method: "jaspar_pwm".into(),
            source: "promoter_scan".into(),
            motif_id: "M1".into(),
            promoter_start: Some(100),
            promoter_end: Some(200),
            hit_start: Some(10),
            hit_end: Some(20),
            strand: "+".into(),
            score: Some(12.5),
            p_value: Some(1e-5),
            ..Default::default()
        }
    }

    #[test]
    fn inserts_exports_and_counts() {
        let path = std::env::temp_dir().join(format!("regutrace_store_{}.db", std::process::id()));
        let _ = std::fs::remove_file(&path);
        {
            let mut store = Store::open(&path).unwrap();
            assert_eq!(store.insert_edges(&[edge()]).unwrap(), 1);
            assert_eq!(store.insert_edges(&[edge()]).unwrap(), 0); // duplicate
            let stats = store.stats().unwrap();
            assert_eq!(stats, vec![("D".into(), "prediction_motif".into(), 1)]);
            let out = store.export_edges().unwrap();
            assert_eq!(out.len(), 1);
            assert_eq!(out[0].target_id, "G1");
        }
        let _ = std::fs::remove_file(path);
    }
}
