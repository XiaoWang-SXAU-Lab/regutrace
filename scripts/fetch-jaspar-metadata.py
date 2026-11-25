#!/usr/bin/env python3
"""Fetch JASPAR matrix metadata for a local PFM file or ID list.

This script only downloads metadata (family, species, UniProt accession,
experimental method, PubMed IDs) from the public JASPAR REST API. It does not
redistribute matrices. JASPAR data are CC BY 4.0; cite the release indicated by
the API when using the output.

Example:
  python scripts/fetch-jaspar-metadata.py \
    --pfm JASPAR2026_CORE_plants_non-redundant_pfms_jaspar.txt \
    --cache ~/.cache/regutrace/jaspar2026 \
    --out jaspar2026-plant-metadata.tsv --workers 8
"""
import argparse
import csv
import json
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

API = "https://jaspar.elixir.no/api/v1/matrix/{matrix_id}/?format=json"


def read_ids(path: Path):
    ids = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            matrix_id = line[1:].split()[0]
            if matrix_id:
                ids.append(matrix_id)
    return ids


def fetch(matrix_id: str, cache: Path, retries: int = 3):
    target = cache / f"{matrix_id}.json"
    if target.is_file():
        return matrix_id, json.loads(target.read_text(encoding="utf-8"))
    last = None
    for attempt in range(retries):
        try:
            request = urllib.request.Request(
                API.format(matrix_id=matrix_id),
                headers={"User-Agent": "regutrace-metadata/1.0"},
            )
            with urllib.request.urlopen(request, timeout=60) as response:
                data = json.load(response)
            target.write_text(json.dumps(data, indent=1), encoding="utf-8")
            return matrix_id, data
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as error:
            last = error
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"{matrix_id}: {last}")


def join(values):
    return ";".join(str(value) for value in values or [] if value not in (None, ""))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--pfm", type=Path)
    source.add_argument("--ids", type=Path)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    ids = read_ids(args.pfm or args.ids)
    if not ids:
        raise SystemExit("no matrix IDs found")
    args.cache.mkdir(parents=True, exist_ok=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    failures = []
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = {pool.submit(fetch, matrix_id, args.cache): matrix_id for matrix_id in ids}
        for index, future in enumerate(as_completed(futures), 1):
            matrix_id = futures[future]
            try:
                _, data = future.result()
                species = data.get("species") or []
                rows.append(
                    {
                        "matrix_id": data.get("matrix_id", matrix_id),
                        "name": data.get("name", ""),
                        "base_id": data.get("base_id", ""),
                        "version": data.get("version", ""),
                        "collection": data.get("collection", ""),
                        "family": join(data.get("family")),
                        "class": join(data.get("class")),
                        "tax_group": data.get("tax_group", ""),
                        "species": join(item.get("name", "") for item in species),
                        "tax_ids": join(item.get("tax_id", "") for item in species),
                        "uniprot_ids": join(data.get("uniprot_ids")),
                        "method": data.get("type", ""),
                        "pubmed_ids": join(data.get("pubmed_ids")),
                        "medline": join(data.get("medline")),
                    }
                )
            except Exception as error:  # noqa: BLE001 - report and continue
                failures.append((matrix_id, str(error)))
            if index % 50 == 0 or index == len(ids):
                print(f"{index}/{len(ids)}", file=sys.stderr)
    rows.sort(key=lambda row: row["matrix_id"])
    columns = [
        "matrix_id",
        "name",
        "base_id",
        "version",
        "collection",
        "family",
        "class",
        "tax_group",
        "species",
        "tax_ids",
        "uniprot_ids",
        "method",
        "pubmed_ids",
        "medline",
    ]
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} metadata rows to {args.out}", file=sys.stderr)
    if failures:
        print(f"{len(failures)} failures:", file=sys.stderr)
        for matrix_id, error in failures[:20]:
            print(f"  {matrix_id}: {error}", file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
