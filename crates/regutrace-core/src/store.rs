//! SQLite evidence store. One row per evidence record; duplicate imports are
//! ignored via a stable record key.

use std::collections::HashSet;
use std::path::Path;

use anyhow::{bail, Context, Result};
use rusqlite::{params, Connection, OptionalExtension};

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
CREATE TABLE IF NOT EXISTS meta (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
"#;

const SCHEMA_VERSION: &str = "2";

pub struct Store {
    conn: Connection,
}

fn opt_u64(v: Option<u64>) -> String {
    v.map(|x| x.to_string()).unwrap_or_default()
}

/// Length-prefix each field so delimiters inside field values cannot make two
/// distinct evidence rows share a key.
fn push_key_part(out: &mut String, value: &str) {
    out.push_str(&value.len().to_string());
    out.push(':');
    out.push_str(value);
    out.push(';');
}

fn checked_coord(name: &str, value: Option<u64>) -> Result<Option<i64>> {
    match value {
        Some(v) => Ok(Some(
            i64::try_from(v).with_context(|| format!("{name} {v} exceeds i64::MAX"))?,
        )),
        None => Ok(None),
    }
}

fn record_key(e: &EdgeRecord) -> String {
    let mut out = String::new();
    for part in [
        &e.tf_id,
        &e.target_id,
        &e.evidence_type,
        &e.evidence_level,
        &e.method,
        &e.source,
        &e.pmid,
        &e.condition,
        &e.motif_id,
    ] {
        push_key_part(&mut out, part);
    }
    for part in [
        opt_u64(e.promoter_start),
        opt_u64(e.promoter_end),
        opt_u64(e.hit_start),
        opt_u64(e.hit_end),
    ] {
        push_key_part(&mut out, &part);
    }
    push_key_part(&mut out, &e.strand);
    out
}

fn migrate_keys(conn: &mut Connection) -> Result<()> {
    let tx = conn.transaction()?;
    let mut rows: Vec<(i64, EdgeRecord)> = Vec::new();
    {
        let mut stmt = tx.prepare(
            "SELECT id, tf_id, tf_symbol, target_id, target_symbol, evidence_type,
                    evidence_level, method, source, pmid, condition, promoter_start,
                    promoter_end, motif_id, hit_start, hit_end, strand, score, p_value, notes
             FROM evidence ORDER BY id",
        )?;
        let iter = stmt.query_map([], |row| {
            Ok((
                row.get::<_, i64>(0)?,
                EdgeRecord {
                    tf_id: row.get(1)?,
                    tf_symbol: row.get(2)?,
                    target_id: row.get(3)?,
                    target_symbol: row.get(4)?,
                    evidence_type: row.get(5)?,
                    evidence_level: row.get(6)?,
                    method: row.get(7)?,
                    source: row.get(8)?,
                    pmid: row.get(9)?,
                    condition: row.get(10)?,
                    promoter_start: row.get::<_, Option<i64>>(11)?.map(|v| v as u64),
                    promoter_end: row.get::<_, Option<i64>>(12)?.map(|v| v as u64),
                    motif_id: row.get(13)?,
                    hit_start: row.get::<_, Option<i64>>(14)?.map(|v| v as u64),
                    hit_end: row.get::<_, Option<i64>>(15)?.map(|v| v as u64),
                    strand: row.get(16)?,
                    score: row.get(17)?,
                    p_value: row.get(18)?,
                    notes: row.get(19)?,
                },
            ))
        })?;
        for r in iter {
            rows.push(r?);
        }
    }
    // Park all keys on unique temporary values first, so assigning a new key
    // cannot collide with a not-yet-migrated row's old key.
    tx.execute("UPDATE evidence SET record_key = 'legacy-' || id", [])?;
    let mut seen: HashSet<String> = HashSet::new();
    for (id, edge) in rows {
        let key = record_key(&edge);
        if seen.insert(key.clone()) {
            tx.execute(
                "UPDATE evidence SET record_key = ?1 WHERE id = ?2",
                params![key, id],
            )?;
        } else {
            tx.execute("DELETE FROM evidence WHERE id = ?1", params![id])?;
        }
    }
    tx.commit()?;
    Ok(())
}

impl Store {
    pub fn open(path: &Path) -> Result<Self> {
        let mut conn = Connection::open(path)
            .with_context(|| format!("open evidence database {}", path.display()))?;
        conn.execute_batch(SCHEMA)?;
        let version: Option<String> = conn
            .query_row(
                "SELECT value FROM meta WHERE key = 'schema_version'",
                [],
                |row| row.get(0),
            )
            .optional()?;
        match version.as_deref() {
            Some(v) if v == SCHEMA_VERSION => {}
            Some(other) => bail!("unsupported evidence database schema version {other}"),
            None => {
                let count: i64 =
                    conn.query_row("SELECT COUNT(*) FROM evidence", [], |row| row.get(0))?;
                if count > 0 {
                    migrate_keys(&mut conn)?;
                }
                conn.execute(
                    "INSERT OR REPLACE INTO meta(key, value) VALUES('schema_version', ?1)",
                    params![SCHEMA_VERSION],
                )?;
            }
        }
        Ok(Self { conn })
    }

    /// Insert evidence rows; returns the number of newly inserted rows.
    pub fn insert_edges(&mut self, edges: &[EdgeRecord]) -> Result<usize> {
        for e in edges {
            if e.tf_id.trim().is_empty() || e.target_id.trim().is_empty() {
                bail!("evidence row missing tf_id or target_id");
            }
            if e.evidence_type.trim().is_empty() {
                bail!("evidence row missing evidence_type");
            }
            for (name, value) in [
                ("promoter_start", e.promoter_start),
                ("promoter_end", e.promoter_end),
                ("hit_start", e.hit_start),
                ("hit_end", e.hit_end),
            ] {
                if value == Some(0) {
                    bail!("{name} must be >= 1");
                }
                let _ = checked_coord(name, value)?;
            }
            if let (Some(start), Some(end)) = (e.promoter_start, e.promoter_end) {
                if start > end {
                    bail!("promoter_start {start} > promoter_end {end}");
                }
            }
            if let (Some(start), Some(end)) = (e.hit_start, e.hit_end) {
                if start > end {
                    bail!("hit_start {start} > hit_end {end}");
                }
            }
            if e.score.is_some_and(|v| !v.is_finite())
                || e.p_value.is_some_and(|v| !v.is_finite())
            {
                bail!(
                    "evidence row has non-finite score or p_value for {} -> {}",
                    e.tf_id,
                    e.target_id
                );
            }
            if let Some(p) = e.p_value {
                if !(0.0..=1.0).contains(&p) {
                    bail!(
                        "p_value {p} out of [0,1] for {} -> {}",
                        e.tf_id,
                        e.target_id
                    );
                }
            }
            if !matches!(e.evidence_level.as_str(), "A" | "B" | "C" | "D" | "E") {
                bail!(
                    "invalid evidence_level '{}' for {} -> {} (expected A/B/C/D/E)",
                    e.evidence_level,
                    e.tf_id,
                    e.target_id
                );
            }
        }
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
                let promoter_start = checked_coord("promoter_start", e.promoter_start)?;
                let promoter_end = checked_coord("promoter_end", e.promoter_end)?;
                let hit_start = checked_coord("hit_start", e.hit_start)?;
                let hit_end = checked_coord("hit_end", e.hit_end)?;
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
                    promoter_start,
                    promoter_end,
                    e.motif_id,
                    hit_start,
                    hit_end,
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

    #[test]
    fn distinct_user_evidence_is_not_collapsed() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_store_user_{}.db",
            std::process::id()
        ));
        let _ = std::fs::remove_file(&path);
        {
            let mut store = Store::open(&path).unwrap();
            let base = EdgeRecord {
                tf_id: "TF1".into(),
                target_id: "G1".into(),
                evidence_type: "point_validation_emsa".into(),
                evidence_level: "A".into(),
                method: "EMSA".into(),
                source: "paper".into(),
                pmid: "111".into(),
                ..Default::default()
            };
            let mut other_pmid = base.clone();
            other_pmid.pmid = "222".into();
            let mut other_level = base.clone();
            other_level.evidence_level = "B".into();
            assert_eq!(
                store
                    .insert_edges(&[base.clone(), other_pmid, other_level])
                    .unwrap(),
                3
            );
            assert_eq!(store.insert_edges(&[base]).unwrap(), 0);
        }
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn migrates_legacy_keys_and_dedups() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_store_migrate_{}.db",
            std::process::id()
        ));
        let _ = std::fs::remove_file(&path);
        let legacy = edge();
        {
            let conn = Connection::open(&path).unwrap();
            conn.execute_batch(SCHEMA).unwrap();
            conn.execute("DELETE FROM meta", []).unwrap();
            let old_key = format!(
                "{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
                legacy.tf_id,
                legacy.target_id,
                legacy.evidence_type,
                legacy.method,
                legacy.source,
                legacy.motif_id,
                legacy.promoter_start.unwrap(),
                legacy.hit_start.unwrap(),
                legacy.strand,
                legacy.condition
            );
            conn.execute(
                "INSERT INTO evidence (
                    record_key, tf_id, tf_symbol, target_id, target_symbol,
                    evidence_type, evidence_level, method, source, pmid, condition,
                    promoter_start, promoter_end, motif_id, hit_start, hit_end, strand,
                    score, p_value, notes
                 ) VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11,?12,?13,?14,?15,?16,?17,?18,?19,?20)",
                params![
                    old_key,
                    legacy.tf_id.as_str(),
                    legacy.tf_symbol.as_str(),
                    legacy.target_id.as_str(),
                    legacy.target_symbol.as_str(),
                    legacy.evidence_type.as_str(),
                    legacy.evidence_level.as_str(),
                    legacy.method.as_str(),
                    legacy.source.as_str(),
                    legacy.pmid.as_str(),
                    legacy.condition.as_str(),
                    legacy.promoter_start.map(|v| v as i64),
                    legacy.promoter_end.map(|v| v as i64),
                    legacy.motif_id.as_str(),
                    legacy.hit_start.map(|v| v as i64),
                    legacy.hit_end.map(|v| v as i64),
                    legacy.strand.as_str(),
                    legacy.score,
                    legacy.p_value,
                    legacy.notes.as_str(),
                ],
            )
            .unwrap();
        }
        let mut store = Store::open(&path).unwrap();
        assert_eq!(store.insert_edges(&[legacy]).unwrap(), 0);
        let version: String = store
            .conn
            .query_row(
                "SELECT value FROM meta WHERE key = 'schema_version'",
                [],
                |row| row.get(0),
            )
            .unwrap();
        assert_eq!(version, SCHEMA_VERSION);
        let count: i64 = store
            .conn
            .query_row("SELECT COUNT(*) FROM evidence", [], |row| row.get(0))
            .unwrap();
        assert_eq!(count, 1);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn rejects_coordinates_above_i64_max() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_store_u64_{}.db",
            std::process::id()
        ));
        let _ = std::fs::remove_file(&path);
        let mut store = Store::open(&path).unwrap();
        let mut bad = edge();
        bad.promoter_end = Some(u64::MAX);
        assert!(store.insert_edges(&[bad]).is_err());
        let count: i64 = store
            .conn
            .query_row("SELECT COUNT(*) FROM evidence", [], |row| row.get(0))
            .unwrap();
        assert_eq!(count, 0);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn rejects_reversed_or_zero_coordinates() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_store_coords_{}.db",
            std::process::id()
        ));
        let _ = std::fs::remove_file(&path);
        let mut store = Store::open(&path).unwrap();
        let mut reversed = edge();
        reversed.promoter_start = Some(20);
        reversed.promoter_end = Some(10);
        assert!(store.insert_edges(&[reversed]).is_err());
        let mut zero = edge();
        zero.hit_start = Some(0);
        assert!(store.insert_edges(&[zero]).is_err());
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn record_key_is_delimiter_safe() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_store_delimiter_{}.db",
            std::process::id()
        ));
        let _ = std::fs::remove_file(&path);
        {
            let mut store = Store::open(&path).unwrap();
            let mut a = edge();
            a.method = "x|y".into();
            a.source = "z".into();
            let mut b = edge();
            b.method = "x".into();
            b.source = "y|z".into();
            assert_eq!(store.insert_edges(&[a, b]).unwrap(), 2);
        }
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn rejects_invalid_level_on_insert() {
        let path = std::env::temp_dir().join(format!(
            "regutrace_store_badlevel_{}.db",
            std::process::id()
        ));
        let _ = std::fs::remove_file(&path);
        let mut store = Store::open(&path).unwrap();
        let mut bad = edge();
        bad.evidence_level = "user_curated".into();
        assert!(store.insert_edges(&[bad]).is_err());
        let _ = std::fs::remove_file(path);
    }
}
