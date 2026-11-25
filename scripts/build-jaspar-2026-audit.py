#!/usr/bin/env python3
"""Prepare a JASPAR donor audit for the plant reference pack.

Step 1 (`fetch`): read the JASPAR metadata TSV, fetch UniProt FASTA for the
accessions attached to one-accession matrices, and write the raw donor table.
Step 2 (`audit`): join the donor table with an `identify-plant-tfs` catalogue
and write the `donor-audit.tsv` consumed by `scripts/build-plant-reference.py`.

This script records provenance; it does not decide that a matrix is a validated
binding model. JASPAR data are CC BY 4.0; UniProt data are distributed under the
UniProt license. Cite both in the pack manifest.
"""
import argparse
import csv
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

UNIPROT = "https://rest.uniprot.org/uniprotkb/search"


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path):
    key = None
    sequence = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if key:
                    yield key, "".join(sequence)
                key = line[1:].strip()
                sequence = []
            elif key:
                sequence.append(line.strip())
    if key:
        yield key, "".join(sequence)


def fetch_uniprot(accessions, cache, workers=5):
    cache.mkdir(parents=True, exist_ok=True)
    def fetch(chunk):
        joined = " OR ".join(f"accession:{accession}" for accession in chunk)
        key = str(abs(hash(joined)))
        target = cache / f"{key}.faa"
        if target.is_file():
            return target.read_text(encoding="utf-8")
        url = UNIPROT + "?" + urllib.parse.urlencode(
            {"query": joined, "format": "fasta", "size": 500}
        )
        last = None
        for attempt in range(4):
            try:
                with urllib.request.urlopen(url, timeout=45) as response:
                    text = response.read().decode()
                target.write_text(text, encoding="utf-8")
                return text
            except (urllib.error.URLError, TimeoutError) as error:
                last = error
                time.sleep(2 * (attempt + 1))
        raise RuntimeError(f"UniProt chunk failed: {last}")
    chunks = [accessions[i : i + 45] for i in range(0, len(accessions), 45)]
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        return "".join(pool.map(fetch, chunks))


def parse_uniprot_header(header):
    parts = header.split("|")
    accession = parts[1] if len(parts) > 1 else header.split()[0]
    taxid = re.search(r"\bOX=(\d+)", header)
    organism = header.split(" OS=", 1)[1].split(" OX=", 1)[0] if " OS=" in header else ""
    return accession, organism, taxid.group(1) if taxid else ""


def step_fetch(args):
    metadata = read_tsv(args.metadata)
    donors = []
    accessions = set()
    for row in metadata:
        ids = [value.strip() for value in row["uniprot_ids"].split(";") if value.strip()]
        if not ids:
            status = "no_accession"
            accession = ""
        elif len(ids) > 1:
            status = "multi_accession_or_complex"
            accession = ""
        else:
            accession = ids[0]
            status = "eligible_pending_classification"
            if re.fullmatch(r"[A-Z0-9]{6,10}", accession):
                accessions.add(accession)
            else:
                status = "non_uniprot_accession"
        donors.append(
            {
                "motif_id": row["matrix_id"],
                "name": row["name"],
                "tax_ids": row["tax_ids"],
                "species": row["species"],
                "accession": accession,
                "method": row["method"],
                "pubmed_ids": row["pubmed_ids"],
                "jaspar_family": row["family"],
                "status": status,
            }
        )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    accessions = sorted(accessions)
    fasta = fetch_uniprot(accessions, args.cache, args.workers)
    original = args.out_dir / "uniprot-original.faa"
    donor_faa = args.out_dir / "jaspar-tfs.faa"
    original.write_text(fasta, encoding="utf-8")
    taxa = {}
    with donor_faa.open("w", encoding="utf-8") as out:
        for header, sequence in read_fasta(original):
            accession, organism, taxid = parse_uniprot_header(header)
            taxa[accession] = (organism, taxid)
            out.write(f">{accession}\n{sequence}\n")
    write_tsv(args.out_dir / "donor-taxa.tsv", [{"accession": key, "organism": value[0], "tax_id": value[1]} for key, value in sorted(taxa.items())], ["accession", "organism", "tax_id"])
    write_tsv(args.out_dir / "jaspar-donors.tsv", donors, list(donors[0]))
    print(f"matrices {len(metadata)}; accessions requested {len(accessions)}; FASTA records {len(taxa)}")
    print(f"wrote {args.out_dir / 'jaspar-donors.tsv'} and {donor_faa}")


def step_audit(args):
    catalogue = {}
    for row in read_tsv(args.catalogue):
        catalogue[row["gene_id"]] = (row["family"], row["identification_status"])
    taxa = {row["accession"]: row for row in read_tsv(args.out_dir / "donor-taxa.tsv")}
    audit = []
    for row in read_tsv(args.out_dir / "jaspar-donors.tsv"):
        accession = row["accession"]
        tax_id = taxa.get(accession, {}).get("tax_id", "")
        family = catalogue.get(accession, ("", ""))[0]
        status = row["status"]
        if status == "eligible_pending_classification":
            if accession not in catalogue:
                status = "accession_not_classified"
            elif catalogue[accession][1] != "identified":
                status = "family_unresolved"
            elif not tax_id or tax_id not in row["tax_ids"].split(";"):
                status = "taxonomy_unresolved"
            elif not family:
                status = "family_unresolved"
            else:
                status = "eligible"
        if row["species"] == "None" and accession in taxa:
            row["species"] = taxa[accession]["organism"]
        row["status"] = status
        row["sequence_tax_id"] = tax_id
        row["planttfhunter_family"] = family
        audit.append(row)
    fields = [
        "motif_id",
        "name",
        "tax_ids",
        "species",
        "accession",
        "status",
        "sequence_tax_id",
        "planttfhunter_family",
        "method",
        "pubmed_ids",
        "jaspar_family",
    ]
    write_tsv(args.out_dir / "donor-audit.tsv", audit, fields)
    from collections import Counter
    print("status", Counter(row["status"] for row in audit))
    print("eligible motifs", sum(row["status"] == "eligible" for row in audit))
    print(f"wrote {args.out_dir / 'donor-audit.tsv'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=["fetch", "audit"])
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--cache", type=Path, default=Path(".jaspar-uniprot-cache"))
    parser.add_argument("--catalogue", type=Path)
    parser.add_argument("--workers", type=int, default=5)
    args = parser.parse_args()
    if args.step == "fetch":
        if not args.metadata:
            parser.error("fetch needs --metadata")
        step_fetch(args)
    else:
        if not args.catalogue:
            parser.error("audit needs --catalogue")
        step_audit(args)


if __name__ == "__main__":
    main()
