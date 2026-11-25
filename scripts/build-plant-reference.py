"""Build a plant reference from a donor audit and a pinned JASPAR PFM library."""
import argparse
import csv
import hashlib
import json
import math
import shutil
from collections import defaultdict
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("legacy", type=Path, help="Plant beta reference containing HMMs and rules")
    parser.add_argument("audit", type=Path, help="Directory with donor-audit.tsv and jaspar-tfs.faa")
    parser.add_argument("output", type=Path, help="New directory for versioned assets")
    parser.add_argument("--motifs", type=Path, help="JASPAR-style PFM file to use instead of the legacy library")
    parser.add_argument("--pack-version", help="manifest pack_version (default: output directory name)")
    parser.add_argument("--strain", help="manifest strain")
    parser.add_argument("--description", help="manifest description")
    parser.add_argument("--motif-scope", help="manifest motif_scope")
    parser.add_argument("--license", help="manifest license")
    parser.add_argument("--attribution", help="manifest attribution")
    args = parser.parse_args()
    manifest = json.loads((args.legacy / "manifest.json").read_text(encoding="utf-8-sig"))
    library = args.motifs or (args.legacy / manifest["files"]["motifs"])
    motif_relative = f"motifs/{library.name}" if args.motifs else manifest["files"]["motifs"]
    information = {}
    for line in library.read_text().splitlines():
        if line.startswith(">"):
            matrix_id = line[1:].split()[0]
            information[matrix_id] = []
        elif "[" in line:
            information[matrix_id].append([float(x) for x in line.split("[", 1)[1].split("]")[0].split()])
    def information_content(matrix_id):
        score = 0.0
        for col in zip(*information[matrix_id]):
            total = sum(col)
            if total:
                score += sum(x / total * math.log2(4 * x / total) for x in col if x)
        return score
    by_accession = defaultdict(list)
    with (args.audit / "donor-audit.tsv").open(newline="", encoding="utf-8-sig") as handle:
        for donor in csv.DictReader(handle, delimiter="\t"):
            if donor["status"] == "eligible":
                by_accession[donor["accession"]].append(donor)
    reps = {
        acc: min(candidates, key=lambda r: (-information_content(r["motif_id"]), r["motif_id"]))
        for acc, candidates in by_accession.items()
    }
    if not reps:
        raise ValueError("No eligible reference donors")
    args.output.mkdir(parents=True, exist_ok=False)
    files = {**manifest["files"], "reference_proteins": "reference-proteins.faa"}
    files["motifs"] = motif_relative
    for key, relative in manifest["files"].items():
        if key in {"reference_map", "motifs"}:
            continue
        target = args.output / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(args.legacy / relative, target)
    motif_target = args.output / motif_relative
    motif_target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(library, motif_target)
    columns = ["tf_id", "tf_symbol", "family", "motif_id", "reference_species", "tax_id"]
    with (args.output / files["reference_map"]).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(columns)
        for acc, r in sorted(reps.items()):
            writer.writerow([acc, r["name"], r["planttfhunter_family"], r["motif_id"], r["species"], r["sequence_tax_id"]])
    recovered = set()
    with (args.audit / "jaspar-tfs.faa").open() as source, (args.output / files["reference_proteins"]).open("w") as dest:
        for line in source:
            if line.startswith(">"):
                acc = line[1:].split()[0]
                keep = acc in reps
                if keep:
                    recovered.add(acc)
            if keep:
                dest.write(line)
    if recovered != reps.keys():
        raise ValueError("Missing reference FASTA records")
    manifest.update(
        pack_version=args.pack_version or args.output.name,
        species="Plants",
        strain=args.strain or "JASPAR 2024 CORE plant donors",
        assembly="Multiple reference species",
        motif_scope=args.motif_scope or "Same-family reciprocal protein best hits; one representative JASPAR matrix per donor",
        description=args.description or "PlantTFHunter family rules and JASPAR 2024 CORE plant motif transfer",
        files=files,
        checksums={name: hashlib.sha256((args.output / name).read_bytes()).hexdigest() for name in files.values()},
    )
    if args.license:
        manifest["license"] = args.license
    if args.attribution:
        manifest["attribution"] = args.attribution
    manifest["generated_by"] = "scripts/build-plant-reference.py"
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (args.output / "SHA256SUMS").write_text("".join(f"{digest}  {name}\n" for name, digest in manifest["checksums"].items()))
    print(f"{len(reps)} donors, {len({r['motif_id'] for r in reps.values()})} representative motifs")


if __name__ == "__main__":
    main()
