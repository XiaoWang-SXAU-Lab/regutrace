# Phase 1 report — reference registry, legacy pack adapter and source audit

Date: 2026-11-24  
Repository: `regutrace`  
Base HEAD: `9ac4837` (plus this working tree; not committed)

## Scope delivered

Phase 1 adds a versioned registry layer without changing mapping or scanning
semantics:

- `regutrace registry inspect` with two input modes:
  - `--reference-pack <dir>` — adapt a legacy pack;
  - `--registry <tsv>` — inspect an existing registry.
- `--motif-library <txt>` — check registry motifs against the packaged PWM
  library and report unused motifs and duplicate matrices;
- `--donor-proteins <faa>` — optional donor FASTA for registry inputs;
- `--summary-json <json>` — per-source, per-family, per-species, QC and library
  summaries.

Registry columns:

```text
motif_id, donor_id, family, source, source_version, accession,
evidence_type, method, species, tax_id, license, donor_sequence_id,
dbd_start, dbd_end, dbd_models, checksum, qc_status
```

`checksum` is `sha256:<hex>` over the uppercase donor residue string returned by
`fasta::read_fasta`; formatting and case are not part of the checksum.

The legacy adapter derives `species` and `tax_id` from optional
`reference_species` / `tax_id` columns in `reference-map.tsv`, falling back to
the manifest species. It does **not** invent licenses, DBD coordinates,
experimental methods or motif provenance that the pack does not contain.

## Validation rules

A pack or registry is rejected when:

- a motif or donor ID is empty;
- the same `(motif_id, donor_id, source, source_version)` appears twice;
- a donor sequence named by the registry is absent from the supplied FASTA;
- a `sha256:` donor checksum does not match the FASTA sequence;
- only one of `dbd_start` / `dbd_end` is present, coordinates are zero/reversed,
  or the end exceeds the donor sequence;
- a manifest checksum does not match the file on disk;
- a file listed in `manifest.files` is missing;
- a registry motif is absent from the supplied PWM library.

`registry inspect --reference-pack` verifies every entry in
`manifest.checksums` before reading the map. Packs without a `checksums` block
are still accepted, but only file existence is checked.

## Real pack result

Input:

```text
plant-jaspar-2024-rbh
manifest pack_version: plant-jaspar-2024-rbh
reference map: 615 donors, one motif each
reference proteins: 615 records
motif library: 805 matrices
```

Command (from the repository root):

```bash
regutrace registry inspect \
  --reference-pack <pack>/plant-jaspar-2024-rbh \
  --motif-library <pack>/plant-jaspar-2024-rbh/motifs/JASPAR2024_CORE_plants.txt \
  --out <out>/plant-jaspar-2024-rbh.registry.tsv \
  --summary-json <out>/plant-jaspar-2024-rbh.registry-summary.json
```

Result:

```text
registry: 615 records, 615 motifs, 615 donors; 615 missing DBD coordinates
motif library: 805 motifs, 190 unused by registry, 0 missing from library,
               3 duplicate normalized-matrix groups
```

Summary:

| item | value |
|---|---|
| records / unique motifs / unique donors | 615 / 615 / 615 |
| families | 34 |
| donor species | 15 |
| QC status | `legacy_unverified` for all 615 |
| license | `unknown` for all 615 |
| DBD coordinates | missing for all 615 |
| donor checksums | present for all 615 |
| library motifs used / unused | 615 / 190 |
| registry motifs missing from library | 0 |
| duplicate normalized-matrix groups | 3 groups / 6 motifs |

Top donor species: Arabidopsis thaliana 541, Zea mays 39, Physcomitrium patens
6, Arabidopsis lyrata 5, Glycine max 5. The complete species and family tables
are in the generated summary JSON.

Exact duplicate matrix groups in the JASPAR 2024 plant CORE library:

| group | motif IDs | donors used in registry | family |
|---|---|---|---|
| 1 | `MA1095.1`, `MA1098.1` | both | TCP |
| 2 | `MA1732.1`, `MA1750.1` | both | E2F/DP |
| 3 | `MA1738.1`, `MA2371.1` | `MA1738.1` only | MYB |

These are byte-identical count matrices, not similarity clusters. Groups 1 and 2
mean two different donors currently contribute the same matrix; Phase 3 must not
count their scan hits as independent TF evidence.

Output hashes:

```text
registry.tsv    f9b3040a411e9544096d3cc35b5092d2b501e9bfcc7ce199e0a667174a917538
summary.json    0caa564db474ea298df7f81f98332a0c58c90b98de49dcf039f774c52f2897fb
```

## Pack provenance audit

### The pack motif library is official JASPAR 2024 CORE plants

The pack file:

```text
motifs/JASPAR2024_CORE_plants.txt
```

has SHA-256:

```text
e25dc968b721cf73ebab4d46bc5432ad80304da7797321bc7637cd319ada0405
```

This is identical to the official file:

```text
https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_plants_non-redundant_pfms_jaspar.txt
```

All 805 matrix IDs in the pack are present in the official file and vice versa
(805 / 805 exact ID overlap). This verifies the motif-library provenance and the
JASPAR 2024 release label in the manifest; it does not verify the donor audit or
the donor-to-motif representative choice.

### JASPAR 2026 is available and licensed CC BY 4.0

Audited on 2026-11-24:

| field | value |
|---|---|
| release | JASPAR 2026, release 11 |
| release date | 2025-12-31 |
| DOI | 10.1093/nar/gkaf1209 |
| PMID | 41325984 |
| license | Creative Commons Attribution 4.0 International (CC BY 4.0) |
| license page | https://jaspar.elixir.no/license/ |
| plant CORE PFM | https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_plants_non-redundant_pfms_jaspar.txt |

The downloaded plant CORE file contains 927 matrices, SHA-256:

```text
617df7f8440fcd92be5527e26485efc7050aca0ddaf5538f6fb909af2d0fd945
```

Difference from the pack's JASPAR 2024 library:

| comparison | count |
|---|---|
| exact shared matrix IDs | 781 |
| JASPAR 2026-only IDs | 146 |
| pack-only (JASPAR 2024) IDs | 24 |
| same base ID but updated version suffix | 5 |

Net change: 805 → 927 matrices. This is an upgrade opportunity, not an
implemented upgrade: the current pack has one representative matrix per donor,
and mapping 615 donor proteins to the expanded 2026 matrix set requires re-running
the donor audit. Swapping the motif file alone would leave 146 new matrices
unmapped and 24 old matrices without a donor.

Downloaded JASPAR files were kept in the local build cache
(`%LOCALAPPDATA%\rsWGCNA-build\jaspar2024` and `jaspar2026`) and were **not**
added to the repository or the sync folder.

### JASPAR 2026 candidate pack built

A complete candidate pack was built from the JASPAR 2026 release and a fresh
donor audit:

```text
pack: plant-jaspar-2026-rbh
manifest SHA-256: c6e2c2aaae1bcb1efdad0ec09323e4fe7cccab3579ac567e410c5ce72fb52e5f
registry SHA-256: bef03c47623e89cc175fdb3bd017a8e464c2a697c6ba957936c41d6151c9b702
pack size: 3.7 MB
```

Build pipeline:

1. Parsed 927 JASPAR 2026 plant CORE matrices from the official PFM file.
2. Fetched metadata for all 927 matrices from the JASPAR REST API; cached
   locally, not committed.
3. Extracted 777 UniProt accessions; fetched 764 FASTA records from UniProt
   release 2026_04 (2026-11-04).
4. Classified donor proteins with the same HMM/rules used by the pipeline:
   751 classified, 733 eligible matrix rows, 713 unique donors.
5. Selected one representative matrix per donor by maximum information content;
   ties broken by motif ID.

Donor audit statuses:

| status | matrices |
|---|---:|
| eligible | 733 |
| no accession | 124 |
| family unresolved | 29 |
| accession not classified | 26 |
| taxonomy unresolved | 11 |
| multiple accessions / complex | 4 |

Pack contents:

| item | value |
|---|---|
| donors / representative motifs | 713 / 713 |
| families | 34 |
| donor species | 15 |
| motif library | 927 JASPAR 2026 plant CORE matrices |
| registry records | 713 |
| registry motifs missing from library | 0 |
| unused library motifs | 214 |
| duplicate normalized-matrix groups | 4 groups / 8 motifs |
| manifest checksums | present and verified |
| embedded registry records | 713 |
| registry license | `CC-BY-4.0` for all 713 |

Registry method breakdown (representative matrices):

```text
DAP-seq 378, PBM 148, JASPAR_CORE 95, ChIP-seq 51, SELEX 26,
ampDAP-seq 9, ChIP-chip 4, SELEX-seq 1, Transfac 1
```

The embedded `registry.tsv` is verified by `registry inspect` against the donor
FASTA checksums. Records are marked `jaspar_metadata_verified`; this means the
JASPAR metadata and UniProt sequence are traceable, not that the donor-to-motif
representative choice has been experimentally validated.

HF1 whole-protein strict mapping with the 2026 pack versus the 2024 pack:

| item | 2024 pack | 2026 pack |
|---|---:|---:|
| mapped TFs | 260 | 271 |
| mapped rows | 260 | 271 |
| common TFs | — | 250 |
| new TFs | — | 21 |
| lost TFs | — | 10 |
| common TFs with changed donor/motif | — | 23 |

Full HF1 pipeline with the 2026 pack (same genome/GFF/CDS/annotation and
promoter settings):

| metric | 2024 pack | 2026 pack |
|---|---:|---:|
| identified TFs | 1,908 | 1,908 |
| mapped TFs | 260 | 271 |
| motifs with hits | 229 | 239 |
| motif hits | 3,068,556 | 3,369,011 |
| evidence rows | 3,090,137 | 3,390,592 |
| TF–target pairs | 1,785,471 | 1,960,433 |
| target union | 37,407 | 37,413 |
| site tie groups | — | 21,581 |
| site tie candidate rows | — | 43,162 |

The 2026 pack increases coverage, but 29 scanned motifs have a minimum
achievable p-value above the default `1e-4` (short matrices such as
`MA2471.1` and `MA2480.1`), so they cannot produce hits at that threshold. This
is a Phase 3 statistical issue, not a pack defect; the scan warns and continues.

The key HF1 site is unchanged: `Fta_6G0039340` still maps to `Q9LTC4` /
`MA0574.1`. Examples of version updates are `MA0566.1 → MA0566.2`,
`MA1378.1 → MA1378.2`, `MA1776.2 → MA1776.3`, and `MA1808.2 → MA1808.3`.

Ten strict mappings are lost. Four are donors whose 2024 matrix was removed
from 2026 (`Q9SDW0/MA1207.2`, `Q9SK03/MA1796.1`, `Q8VZ20/MA1733.2`,
`Q9LVG2/MA2363.1`). Six donors remain in the 2026 pack but no longer win the
reciprocal best-hit competition against the expanded donor set
(`Q39204/MA0566.2`, `O49403/MA1759.2`, `Q9SSA8/MA1256.2`,
`Q9LW31/MA0565.3`, `Q9LUH8/MA1664.3`, `Q9SI37/MA0589.2`). These differences
need target-level benchmarking before the 2026 pack can become the default.

The pack is stored at:

```text
analysis/regutrace-yeast-validation/reference-pack/plant-jaspar-2026-rbh/
```

It includes `NOTICE.txt` with the JASPAR 2026 CC BY 4.0 attribution, the
UniProt 2026_04 CC BY 4.0 attribution, and the build method. It is a candidate
pack, not the default. No external matrix data was committed to the source
repository.

### CIS-BP is accessible but not practical or licensed for redistribution

CIS-BP 3.10 direct URLs exist and respond to a HEAD request:

```text
https://cisbp.ccbr.utoronto.ca/data/3_10/DataFiles/Bulk_downloads/EntireDataset/TF_Information.zip
https://cisbp.ccbr.utoronto.ca/data/3_10/DataFiles/Bulk_downloads/EntireDataset/PWMs.zip
```

`TF_Information.zip` reports `Content-Length: 240452250` and
`Last-Modified: Thu, 28 May 2026`. From this network, GET transfers were
throttled to roughly 7.5 KB in 120 s and then timed out; the complete archives
were not downloaded. The official bulk page and citation page were also slow or
reset.

The RSAT mirror is reachable and small enough to inspect:

```text
https://rsat.eead.csic.es/plants/motif_databases/cisBP/
https://rsat.eead.csic.es/plants/motif_databases/cisBP/CisBp_A_thaliana.tf
```

`CisBp_A_thaliana.tf` is 231,312 bytes in TRANSFAC-like format. Its
`LICENSE.txt` only links PubMed 25215497 and makes a general statement about
open-access articles; it does **not** grant redistribution rights. re3data
labels the repository "Public Domain", but that is metadata about the registry
record, not a content license.

Status: **accessible for local analysis; redistribution blocked**. Do not add
CIS-BP matrices to a distributed pack until written terms are obtained.

### PlantTFDB is accessible; license not stated

The experimental-motif information table and archive are readable:

```text
https://planttfdb.gao-lab.org/download/motif/PlantTFDB_TF_binding_motifs_from_experiments_information.txt
https://planttfdb.gao-lab.org/download/motif/PlantTFDB_TF_binding_motifs_from_experiments.gz
```

The information table is 45,380 bytes / 675 lines. The gzip archive is
118,129 bytes and contains 674 MEME-format motifs (`MOTIF MP00001 ...`). The
site has a "How to Cite" section but no explicit data license or redistribution
statement was found.

Status: **accessible for local analysis; redistribution unverified**. An
adapter can be built, but the matrices should not be redistributed until the
terms are confirmed.

### Other sources

PlantRegMap is partially reachable but not audited. Plant Cistrome and PlantPAN
5.0 were not audited in this phase; PlantPAN's standalone Docker package is
reported to contain 2,254 old motifs rather than the full 7,762 and lacks the
complete gene-ID/protein annotation, so it is a low-priority supplement. No data
from these sources was downloaded.

## DBD coordinate finding

The pack contains `reference-domains.faa` with 138 records, all named as
Arabidopsis AGI IDs (`AT1G06180`, `AT1G08810`, ...). None of these IDs matches
the 615 donor accessions in `reference-proteins.faa` (0 overlap). The file is a
family/DBD model reference, not per-donor DBD annotation.

Consequences:

- the registry cannot fill `dbd_start` / `dbd_end` from this file;
- `dbd_models.json` contains only `{"MYB": ["PF00249", "PF00249"]}`;
- Phase 2 must extract donor DBDs from `reference-proteins.faa` with the same
  HMMER models and rules used for query TFs, and record the coordinate system
  explicitly.

This is the main technical blocker between Phase 1 and DBD-based transfer.

## Tests run

Environment on Windows:

```bash
export PATH="$LOCALAPPDATA/w64devkit/w64devkit/bin:$PATH"
export CC="$LOCALAPPDATA/w64devkit/w64devkit/bin/gcc.exe"
export CXX="$LOCALAPPDATA/w64devkit/w64devkit/bin/g++.exe"
export RUSTFLAGS="-C link-self-contained=yes"
cargo +stable-x86_64-pc-windows-gnu test --release --locked
```

Result:

```text
36 CLI unit tests passed (Phase 1 snapshot)
5 CLI integration tests passed (Phase 1 snapshot)
69 core unit tests passed
0 failed
```

Current suite after the Phase 2/2.5 additions: 40 CLI unit + 7 CLI integration
+ 69 core, all passing.

Phase 1 tests:

- `registry_reads_and_validates_provenance` — TSV schema, DBD pairing, summary;
- `registry_rejects_unpaired_or_invalid_dbd_coordinates`;
- `registry_rejects_duplicate_records_and_missing_donor_sequences`;
- `legacy_pack_adapter_reads_per_donor_species` — species, tax_id, checksum,
  manifest checksum mismatch;
- `motif_library_reports_unused_and_duplicate_matrices` — exact duplicates,
  unused motifs, hard failure for a motif absent from the library;
- `registry_inspect_accepts_legacy_pack_and_writes_summary` — end-to-end CLI,
  including the `inspect` subcommand and unknown-subcommand rejection;
- core PWM tests cover malformed matrices and duplicate motif IDs
  (`pwm::tests::rejects_malformed_motif_tokens`,
  `pwm::tests::rejects_duplicate_motif_ids`, and related normalization tests).

No HF1 rerun was performed in Phase 1; no mapping or scan output changed.

## Compatibility

- `registry` is a new command; existing commands and outputs are unchanged.
- The registry TSV is additive to the pack contract; packs do not yet consume it.
- Manifest checksum verification is strict when a `checksums` block exists.
  All current packs must therefore keep their checksums in sync with files.
- Legacy packs without per-donor DBD coordinates are accepted and marked
  `legacy_unverified`.
- No database migration is required.

## Phase 1 acceptance status

| criterion | status |
|---|---|
| versioned registry schema with provenance fields | met (legacy scope) |
| legacy pack adapter and validation | met |
| per-source / per-family coverage | met |
| duplicate statistics and unused-motif exclusion | met |
| old-pack compatibility | met |
| pack checksum verification | met |
| pack rebuild from registry | **not met** — audit inputs absent |
| JASPAR 2026 integration | **blocked** — license OK, donor re-audit absent |
| CIS-BP integration | **blocked** — license terms unavailable |
| other sources | **blocked** — not audited |

Phase 1 is complete for the legacy registry scope and explicitly partial for the
external-source scope. The pack is verifiable and offline-usable, but a full
rebuild cannot be claimed from this repository: `scripts/build-plant-reference.py`
requires `donor-audit.tsv` and `jaspar-tfs.faa`, which are not present in the
repository.

## Next phase

Phase 2 should start with donor DBD extraction, not with new motif sources:

1. run the existing HMMER/rules pipeline on `reference-proteins.faa`;
2. store per-donor DBD coordinates, models, architecture and completeness;
3. implement DBD-based transfer behind a new mode, leaving legacy strict
   behaviour unchanged;
4. version all thresholds and record their source; abstain where no defensible
   threshold exists.
