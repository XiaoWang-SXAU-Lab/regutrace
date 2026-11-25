#!/usr/bin/env python3
"""Build the Phase 2.5 calibration donor/motif dataset.

Inputs are the existing JASPAR 2024 and JASPAR 2026 plant reference assets and
the DBD extractions produced by `regutrace dbd-extract`. The script merges
releases by donor and motif, records provenance, removes exact duplicate
donor-motif records, groups exact PWM duplicates, and writes exclusions.
"""
import argparse
import collections
import csv
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    information_content,
    matrix_hash,
    method_is_experimental,
    read_2024_metadata,
    read_2026_metadata,
    read_fasta,
    read_pfm,
    read_tsv,
    write_tsv,
)

MIN_DONORS = 15
MIN_SPECIES = 3
MIN_EXACT_GROUPS = 2


def dbd_rows(path):
    rows = {}
    for row in read_tsv(path):
        rows[row["sequence_id"]] = row
    return rows


def dbd_fasta(path):
    return {name: sequence for name, sequence in read_fasta(path)}


def reference_maps(*paths):
    donors = {}
    for path in paths:
        for row in read_tsv(path):
            donor = row["tf_id"]
            entry = donors.setdefault(
                donor,
                {
                    "donor_id": donor,
                    "family": row["family"],
                    "species": row.get("reference_species", ""),
                    "tax_id": row.get("tax_id", ""),
                    "sources": set(),
                },
            )
            if entry["family"] != row["family"]:
                entry["family_conflict"] = f"{entry['family']}|{row['family']}"
            if entry["species"] and row.get("reference_species") and entry["species"] != row["reference_species"]:
                entry["species_conflict"] = f"{entry['species']}|{row['reference_species']}"
            entry["sources"].add(str(path.parent.name))
    return donors


def audit_motifs(path, metadata, release):
    records = []
    exclusions = []
    for row in read_tsv(path):
        donor = row.get("accession", "").strip()
        motif = row.get("motif_id", "").strip()
        if not donor or not motif:
            continue
        if row.get("status") != "eligible":
            exclusions.append(
                {
                    "donor_id": donor,
                    "motif_id": motif,
                    "family": row.get("planttfhunter_family", ""),
                    "reason": f"{release}_audit_{row.get('status', 'unknown')}",
                    "detail": "",
                }
            )
            continue
        meta = metadata.get(motif)
        method = (meta or {}).get("method", row.get("method", ""))
        if not method_is_experimental(method):
            exclusions.append(
                {
                    "donor_id": donor,
                    "motif_id": motif,
                    "family": row.get("planttfhunter_family", ""),
                    "reason": f"{release}_motif_not_experimental",
                    "detail": method or "empty_method",
                }
            )
            continue
        records.append(
            {
                "donor_id": donor,
                "motif_id": motif,
                "family": row.get("planttfhunter_family", ""),
                "release": release,
                "method": method,
                "species": row.get("species", ""),
                "tax_id": row.get("sequence_tax_id", ""),
            }
        )
    return records, exclusions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pack2024", type=Path, required=True)
    parser.add_argument("--metadata2024db", type=Path, required=True)
    parser.add_argument("--audit2024", type=Path, required=True)
    parser.add_argument("--dbd2024tsv", type=Path, required=True)
    parser.add_argument("--dbd2024faa", type=Path, required=True)
    parser.add_argument("--pack2026", type=Path, required=True)
    parser.add_argument("--metadata2026", type=Path, required=True)
    parser.add_argument("--audit2026", type=Path, required=True)
    parser.add_argument("--dbd2026tsv", type=Path, required=True)
    parser.add_argument("--dbd2026faa", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    pfm2024 = read_pfm(args.pack2024 / "motifs/JASPAR2024_CORE_plants.txt")
    pfm2026 = read_pfm(args.pack2026 / "motifs/JASPAR2026_CORE_plants_non-redundant_pfms_jaspar.txt")
    metadata2024 = read_2024_metadata(args.metadata2024db, pfm2024.keys())
    metadata2026 = read_2026_metadata(args.metadata2026)

    # Matrix provenance: union releases; shared IDs must have identical matrices.
    matrix_by_motif = {}
    motif_sources = collections.defaultdict(set)
    matrix_conflicts = []
    for release, matrices in (("JASPAR2024", pfm2024), ("JASPAR2026", pfm2026)):
        for motif, matrix in matrices.items():
            if motif in matrix_by_motif and matrix_by_motif[motif] != matrix:
                matrix_conflicts.append(motif)
            matrix_by_motif[motif] = matrix
            motif_sources[motif].add(release)

    hash_groups = collections.defaultdict(list)
    for motif, matrix in matrix_by_motif.items():
        hash_groups[matrix_hash(matrix)].append(motif)
    exact_group = {}
    exact_group_rows = []
    for digest, motifs in sorted(hash_groups.items()):
        motifs.sort()
        group_id = f"EX{len(exact_group_rows) + 1:04d}"
        for motif in motifs:
            exact_group[motif] = group_id
        exact_group_rows.append(
            {
                "exact_group": group_id,
                "pwm_hash": digest,
                "motifs": ";".join(motifs),
                "size": len(motifs),
            }
        )
    write_tsv(args.out / "exact-pwm-groups.tsv", exact_group_rows)

    donor_records = {}
    exclusions = []
    for release, audit, metadata in (
        ("JASPAR2024", args.audit2024, metadata2024),
        ("JASPAR2026", args.audit2026, metadata2026),
    ):
        records, rejected = audit_motifs(audit, metadata, release)
        exclusions.extend(rejected)
        for record in records:
            donor = donor_records.setdefault(
                record["donor_id"],
                {
                    "donor_id": record["donor_id"],
                    "family": record["family"],
                    "motifs": collections.defaultdict(set),
                    "species": set(),
                    "tax_ids": set(),
                    "releases": set(),
                    "methods": collections.defaultdict(set),
                },
            )
            if donor["family"] != record["family"]:
                donor["family_conflict"] = f"{donor['family']}|{record['family']}"
            donor["motifs"][record["motif_id"]].add(record["release"])
            if record["species"]:
                donor["species"].add(record["species"])
            if record["tax_id"]:
                donor["tax_ids"].add(record["tax_id"])
            donor["releases"].add(release)
            donor["methods"][record["motif_id"]].add(record["method"])

    references = reference_maps(
        args.pack2024 / "reference-map.tsv", args.pack2026 / "reference-map.tsv"
    )
    dbd2024 = dbd_rows(args.dbd2024tsv)
    dbd2026 = dbd_rows(args.dbd2026tsv)
    fasta2024 = dbd_fasta(args.dbd2024faa)
    fasta2026 = dbd_fasta(args.dbd2026faa)

    donor_rows = []
    for donor_id, donor in sorted(donor_records.items()):
        ref = references.get(donor_id, {})
        family = donor["family"] or ref.get("family", "")
        if donor.get("family_conflict") or ref.get("family_conflict"):
            exclusions.append(
                {
                    "donor_id": donor_id,
                    "motif_id": "",
                    "family": family,
                    "reason": "family_conflict",
                    "detail": donor.get("family_conflict") or ref.get("family_conflict", ""),
                }
            )
            continue
        if ref.get("species_conflict"):
            exclusions.append(
                {
                    "donor_id": donor_id,
                    "motif_id": "",
                    "family": family,
                    "reason": "species_conflict",
                    "detail": ref["species_conflict"],
                }
            )
            continue
        if donor_id in dbd2026:
            dbd = dbd2026[donor_id]
            sequence = fasta2026.get(donor_id)
            dbd_source = "JASPAR2026"
        elif donor_id in dbd2024:
            dbd = dbd2024[donor_id]
            sequence = fasta2024.get(donor_id)
            dbd_source = "JASPAR2024"
        else:
            exclusions.append(
                {
                    "donor_id": donor_id,
                    "motif_id": "",
                    "family": family,
                    "reason": "no_dbd_sequence",
                    "detail": "family absent or unsupported by DBD model config",
                }
            )
            continue
        if dbd.get("mode") != "full_architecture":
            exclusions.append(
                {
                    "donor_id": donor_id,
                    "motif_id": "",
                    "family": family,
                    "reason": "dbd_not_full_architecture",
                    "detail": dbd.get("mode", ""),
                }
            )
            continue
        if not sequence:
            exclusions.append(
                {
                    "donor_id": donor_id,
                    "motif_id": "",
                    "family": family,
                    "reason": "missing_dbd_fasta",
                    "detail": "",
                }
            )
            continue
        motifs = sorted(donor["motifs"])
        motif_methods = []
        for motif in motifs:
            motif_methods.extend(sorted(donor["methods"][motif]))
        methods = sorted(set(motif_methods))
        donor_rows.append(
            {
                "donor_id": donor_id,
                "family": family,
                "species": ";".join(sorted(donor["species"] or {ref.get("species", "")} - {""})),
                "tax_id": ";".join(sorted(donor["tax_ids"] or {ref.get("tax_id", "")} - {""})),
                "dbd_start": dbd.get("dbd_start", ""),
                "dbd_end": dbd.get("dbd_end", ""),
                "domain_count": dbd.get("domain_count", ""),
                "dbd_mode": dbd.get("mode", ""),
                "dbd_source": dbd_source,
                "dbd_length": len(sequence),
                "dbd_sequence": sequence,
                "n_motifs": len(motifs),
                "motif_ids": ";".join(motifs),
                "methods": ";".join(methods),
                "source_releases": ";".join(sorted(donor["releases"])),
            }
        )

    motif_records = collections.defaultdict(lambda: {"donors": set(), "methods": set()})
    for donor in donor_rows:
        for motif in donor["motif_ids"].split(";"):
            if not motif:
                continue
            motif_records[motif]["donors"].add(donor["donor_id"])
            for method in donor["methods"].split(";"):
                if method:
                    motif_records[motif]["methods"].add(method)

    motif_rows = []
    for motif, record in sorted(motif_records.items()):
        matrix = matrix_by_motif[motif]
        meta2026 = metadata2026.get(motif, {})
        meta2024 = metadata2024.get(motif, {})
        method = meta2026.get("method") or meta2024.get("method") or ""
        family = meta2026.get("family") or meta2024.get("family") or ""
        species = meta2026.get("species") or meta2024.get("species") or []
        uniprot = meta2026.get("uniprot_ids") or meta2024.get("uniprot_ids") or []
        pubmed = meta2026.get("pubmed_ids") or meta2024.get("pubmed_ids") or []
        motif_rows.append(
            {
                "motif_id": motif,
                "pwm_hash": matrix_hash(matrix),
                "exact_group": exact_group.get(motif, ""),
                "length": len(matrix),
                "information_content": f"{information_content(matrix):.4f}",
                "jaspar_family": family,
                "source_releases": ";".join(sorted(motif_sources[motif])),
                "methods": ";".join(sorted(record["methods"])),
                "method_metadata": method,
                "pubmed_ids": ";".join(pubmed),
                "species": ";".join(item.get("name", "") for item in species),
                "tax_ids": ";".join(str(item.get("tax_id", "")) for item in species),
                "uniprot_ids": ";".join(uniprot),
                "donor_count": len(record["donors"]),
                "donor_ids": ";".join(sorted(record["donors"])),
            }
        )

    family = collections.defaultdict(lambda: {"donors": [], "motifs": set(), "species": set()})
    for donor in donor_rows:
        entry = family[donor["family"]]
        entry["donors"].append(donor)
        entry["motifs"].update(donor["motif_ids"].split(";"))
        entry["species"].update(donor["species"].split(";"))
    family_rows = []
    for name, entry in sorted(family.items()):
        donors = entry["donors"]
        exact_groups = {exact_group.get(motif, motif) for motif in entry["motifs"] if motif}
        n_motifs = len(entry["motifs"] - {""})
        status = (
            "supported"
            if len(donors) >= MIN_DONORS
            and len(entry["species"] - {""}) >= MIN_SPECIES
            and len(exact_groups) >= MIN_EXACT_GROUPS
            else "abstain"
        )
        reason = ""
        if status == "abstain":
            reasons = []
            if len(donors) < MIN_DONORS:
                reasons.append(f"donors<{MIN_DONORS}")
            if len(entry["species"] - {""}) < MIN_SPECIES:
                reasons.append(f"species<{MIN_SPECIES}")
            if len(exact_groups) < MIN_EXACT_GROUPS:
                reasons.append(f"exact_groups<{MIN_EXACT_GROUPS}")
            reason = ";".join(reasons)
        family_rows.append(
            {
                "family": name,
                "donors": len(donors),
                "species": len(entry["species"] - {""}),
                "motifs": n_motifs,
                "exact_groups": len(exact_groups),
                "median_motifs_per_donor": f"{statistics.median(d['n_motifs'] for d in donors):.1f}",
                "median_dbd_length": f"{statistics.median(d['dbd_length'] for d in donors):.1f}",
                "calibration_status": status,
                "abstain_reason": reason,
            }
        )

    donor_rows.sort(key=lambda row: (row["family"], row["donor_id"]))
    write_tsv(
        args.out / "calibration-donors.tsv",
        [{key: value for key, value in row.items() if key != "dbd_sequence"} for row in donor_rows],
    )
    write_tsv(args.out / "calibration-motifs.tsv", motif_rows)
    write_tsv(args.out / "calibration-family-summary.tsv", family_rows)
    write_tsv(args.out / "calibration-exclusions.tsv", exclusions)
    with (args.out / "calibration-dbd.fasta").open("w", encoding="utf-8") as handle:
        for row in donor_rows:
            handle.write(f">{row['donor_id']}\n{row['dbd_sequence']}\n")

    protein_sequences = {}
    for path in (
        args.pack2024 / "reference-proteins.faa",
        args.pack2026 / "reference-proteins.faa",
    ):
        if path.is_file():
            for donor_id, sequence in read_fasta(path):
                protein_sequences.setdefault(donor_id, sequence)
    with (args.out / "calibration-proteins.faa").open("w", encoding="utf-8") as handle:
        for row in donor_rows:
            sequence = protein_sequences.get(row["donor_id"])
            if sequence is None:
                raise SystemExit(f"missing full-length protein for {row['donor_id']}")
            handle.write(f">{row['donor_id']}\n{sequence}\n")

    print(f"calibration donors: {len(donor_rows)}")
    print(f"calibration motifs: {len(motif_rows)}")
    print(f"exact PWM groups: {len(exact_group_rows)}")
    print(f"exclusions: {len(exclusions)}")
    print(f"matrix conflicts: {len(matrix_conflicts)}")
    supported = [row["family"] for row in family_rows if row["calibration_status"] == "supported"]
    print(f"supported families: {len(supported)} / {len(family_rows)}")
    print("supported:", ", ".join(sorted(supported)))


if __name__ == "__main__":
    main()
