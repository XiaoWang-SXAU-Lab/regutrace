"""Shared parsing and PWM helpers for Phase 2.5 calibration."""
import csv
import hashlib
import json
import math
import re
import sqlite3
from pathlib import Path

BASES = "ACGT"
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
EXPERIMENTAL_METHODS = {
    "DAP-seq",
    "ampDAP-seq",
    "PBM",
    "ChIP-seq",
    "ChIP-chip",
    "ChIP-exo",
    "SELEX",
    "SELEX-seq",
    "HT-SELEX",
    "CAP-SELEX",
    "bacterial 1-hybrid",
    "PBM, CSA and/or DIP-chip",
}


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path):
    name = None
    sequence = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name:
                    yield name, "".join(sequence)
                name = line[1:].split()[0]
                sequence = []
            elif name:
                sequence.append(line.strip())
    if name:
        yield name, "".join(sequence)


def read_pfm(path):
    """Read a JASPAR-style PFM file and normalize each column to probabilities.

    Uses the same 0.5 pseudocount policy as the Rust PWM reader for count
    matrices, and column-sum normalization for probability matrices.
    """
    motifs = {}
    motif_id = None
    rows = []
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.strip()
        if line.startswith(">"):
            if motif_id:
                motifs[motif_id] = normalize_rows(rows)
            motif_id = line[1:].split()[0]
            rows = []
        elif motif_id and "[" in line and "]" in line:
            values = line.split("[", 1)[1].split("]", 1)[0].split()
            rows.append([float(value) for value in values])
    if motif_id:
        motifs[motif_id] = normalize_rows(rows)
    return motifs


def normalize_rows(rows):
    if len(rows) != 4:
        raise ValueError(f"expected four base rows, found {len(rows)}")
    length = len(rows[0])
    if any(len(row) != length for row in rows):
        raise ValueError("unequal PWM row lengths")
    columns = []
    for index in range(length):
        column = [rows[base][index] for base in range(4)]
        total = sum(column)
        if max(column) <= 1.0 + 1e-9 and abs(total - 1.0) < 0.15:
            denominator = total + 4e-3
            normalized = [(value + 1e-3) / denominator for value in column]
        else:
            denominator = total + 2.0
            normalized = [(value + 0.5) / denominator for value in column]
        columns.append(normalized)
    return columns


def reverse_complement_matrix(matrix):
    out = []
    for column in reversed(matrix):
        out.append([column[3], column[2], column[1], column[0]])
    return out


def matrix_hash(matrix, digits=6):
    payload = "|".join(
        ",".join(f"{value:.{digits}f}" for value in column) for column in matrix
    )
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def information_content(matrix):
    score = 0.0
    for column in matrix:
        for value in column:
            if value > 0:
                score += value * math.log2(4 * value)
    return score


def pearson(matrix_a, matrix_b, offset):
    """Pearson correlation over overlapping columns at a given offset.

    Offset > 0 shifts B right relative to A; offset < 0 shifts B left.
    Returns (correlation, overlap) or (None, 0).
    """
    a_start = max(0, -offset)
    b_start = max(0, offset)
    overlap = min(len(matrix_a) - a_start, len(matrix_b) - b_start)
    if overlap < 3:
        return None, overlap
    xs = []
    ys = []
    for index in range(overlap):
        xs.extend(matrix_a[a_start + index])
        ys.extend(matrix_b[b_start + index])
    n = len(xs)
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    den_x = sum((x - mean_x) ** 2 for x in xs)
    den_y = sum((y - mean_y) ** 2 for y in ys)
    if den_x <= 1e-12 or den_y <= 1e-12:
        return None, overlap
    return num / math.sqrt(den_x * den_y), overlap


def motif_similarity(matrix_a, matrix_b):
    """TOMTOM-style maximum Pearson correlation over offsets and strands.

    This is the correlation component of TOMTOM, not its E-value. It is used as
    an operational similarity measure, never as biological equivalence proof.
    """
    best = {"correlation": -1.0, "offset": 0, "overlap": 0, "orientation": "forward"}
    for orientation, matrix in (("forward", matrix_b), ("reverse", reverse_complement_matrix(matrix_b))):
        for offset in range(-len(matrix_a) + 3, len(matrix) - 2):
            correlation, overlap = pearson(matrix_a, matrix, offset)
            if correlation is None:
                continue
            if correlation > best["correlation"]:
                best = {
                    "correlation": correlation,
                    "offset": offset,
                    "overlap": overlap,
                    "orientation": orientation,
                }
    if best["correlation"] < 0:
        return None
    return best


def same_genus(species_a, species_b):
    return species_a.split()[0].lower() == species_b.split()[0].lower()


def read_2024_metadata(db_path, motif_ids):
    """Read JASPAR 2024 CORE matrix metadata from the local SQLite dump."""
    connection = sqlite3.connect(db_path)
    placeholders = ",".join("?" for _ in motif_ids)
    rows = connection.execute(
        f"""
        SELECT m.ID, m.BASE_ID || '.' || m.VERSION AS motif_id, m.NAME,
               MAX(CASE WHEN a.TAG='type' THEN a.VAL END) AS method,
               MAX(CASE WHEN a.TAG='family' THEN a.VAL END) AS family,
               MAX(CASE WHEN a.TAG='class' THEN a.VAL END) AS class,
               MAX(CASE WHEN a.TAG='source' THEN a.VAL END) AS source,
               MAX(CASE WHEN a.TAG='tax_group' THEN a.VAL END) AS tax_group
        FROM MATRIX m
        LEFT JOIN MATRIX_ANNOTATION a ON a.ID = m.ID
        WHERE m.COLLECTION='CORE'
          AND (m.BASE_ID || '.' || m.VERSION) IN ({placeholders})
        GROUP BY m.ID
        """,
        list(motif_ids),
    ).fetchall()
    by_motif = {}
    for row in rows:
        matrix_id, motif_id, name, method, family, klass, source, tax_group = row
        proteins = [
            value
            for (value,) in connection.execute(
                "SELECT ACC FROM MATRIX_PROTEIN WHERE ID=?", (matrix_id,)
            )
            if value
        ]
        species = []
        for (tax_id,) in connection.execute(
            "SELECT TAX_ID FROM MATRIX_SPECIES WHERE ID=?", (matrix_id,)
        ):
            name_row = connection.execute(
                "SELECT SPECIES FROM TAX WHERE TAX_ID=?", (tax_id,)
            ).fetchone()
            species.append(
                {"tax_id": str(tax_id), "name": name_row[0] if name_row else ""}
            )
        by_motif[motif_id] = {
            "matrix_id": motif_id,
            "name": name,
            "method": method or "",
            "family": family or "",
            "class": klass or "",
            "source": source or "",
            "tax_group": tax_group or "",
            "uniprot_ids": proteins,
            "species": species,
            "pubmed_ids": [],
        }
    connection.close()
    return by_motif


def read_2026_metadata(path):
    rows = read_tsv(path)
    return {
        row["matrix_id"]: {
            "matrix_id": row["matrix_id"],
            "name": row["name"],
            "method": row["method"],
            "family": row["family"],
            "class": row["class"],
            "source": "JASPAR",
            "tax_group": row["tax_group"],
            "uniprot_ids": [value for value in row["uniprot_ids"].split(";") if value],
            "species": [
                {"tax_id": tax_id, "name": species}
                for tax_id, species in zip(
                    row["tax_ids"].split(";"), row["species"].split(";")
                )
            ],
            "pubmed_ids": [value for value in row["pubmed_ids"].split(";") if value],
        }
        for row in rows
    }


def method_is_experimental(method):
    if not method:
        return False
    return method.strip() in EXPERIMENTAL_METHODS
