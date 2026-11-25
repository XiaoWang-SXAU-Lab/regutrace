# Phase 2 report — DBD extraction and exploratory DBD-based transfer

Date: 2026-11-24  
Repository: `regutrace`  
Base HEAD: `9ac4837` (plus this working tree; not committed)

## Scope delivered

Phase 2 adds two commands and a versioned DBD-model configuration, without
changing the legacy strict mapping path:

- `dbd-extract` — extract one DNA-binding-domain sequence per TF gene from a
  catalogue plus protein FASTA;
- `map-dbd-motifs` — align query DBDs against donor DBDs and emit a
  TF–motif candidate table with explicit provenance and warning labels;
- `config/dbd-models/plant-v1.json` — family-to-DBD-model architecture config;
- `config/dbd-thresholds/exploratory-v1.json` — explicitly uncalibrated
  engineering thresholds used only for pipeline validation.

New registry/configuration fields are additive. `map-tf-motifs`, `tf_map.tsv`
and scan outputs are unchanged.

## DBD model policy

A DBD is extracted only when a family has a single, documented architecture:

- one model or a fixed repeated model (`MYB` = two `PF00249`);
- a fixed ordered pair (`RAV` = `PF00847` + `PF02362`);
- a single domain when the pack rule names an auxiliary domain that is excluded
  (`ARF` uses `PF02362`, not `PF06507`; `GRF` uses `PF08880`, not `PF08879`;
  `HB-PHD` uses `PF00046`, not `PF00628`).

Families with variable arrays or composite DNA binding are marked `abstain` and
produce no DBD sequence:

```text
C2H2, C3H, CO-like, CPP, GRAS, LSD, NF-YA, NF-YB, NF-YC, ZF-HD
```

This is deliberate: the extractor does not invent a boundary for a zinc-finger
array or for a heterotrimeric NF-Y complex. C2H2/C3H array extraction is future
work, not a silent fallback.

The model source is the Pfam/InterPro domain identity used by the plant rule
pack. The config is versioned (`plant-dbd-models-v1`) and each supported family
records its source and excluded auxiliary domains.

## Real-pack donor extraction

Input:

```text
reference proteins: plant-jaspar-2024-rbh/reference-proteins.faa (615)
HMM:                hmm/combined-cli.hmm (64 models)
rules:              rules/planttfdb-reference.json
config:             config/dbd-models/plant-v1.json
```

Step 1 — classify donor proteins with the pack:

```bash
regutrace identify-plant-tfs \
  --proteins reference-proteins.faa \
  --hmm combined-cli.hmm \
  --rules planttfdb-reference.json \
  --gene-ids --hmmsearch <hmmsearch> --threads 4 \
  --out donor-tf_genes.tsv
```

Result: 615/615 donor proteins classified; all `identified`. Family labels
matched the reference map exactly (0 mismatches).

Step 2 — extract donor DBDs:

```bash
regutrace dbd-extract \
  --tf-genes donor-tf_genes.tsv \
  --proteins reference-proteins.faa \
  --dbd-models config/dbd-models/plant-v1.json \
  --out-faa donor-dbd.faa \
  --out-tsv donor-dbd.tsv
```

Result:

| item | value |
|---|---|
| donor proteins | 615 |
| DBD sequences extracted | 585 |
| abstained donors | 30 |
| supported families in config | 45 |
| abstained families in config | 10 |
| extraction mode | `full_architecture` for all 585 |
| coordinate/sequence mismatches | 0 |

The 30 missing donors are exactly the abstained families present in this pack:
C2H2 (21), CPP (4), ZF-HD (5). Every extracted sequence matches the full-protein
slice at its recorded coordinates; every donor DBD row has a `sha256:` checksum.

## Exploratory HF1 run

HF1 inputs used only for engineering comparison:

```text
catalogue: hf1-plant-rbh-acceptance-2026-11-17/pep/tf_genes.tsv
proteins:  _rsWGCNA-local/hf1-orthology-motif-20261117/hf1-tfs.faa
strict map:hf1-plant-rbh-acceptance-2026-11-17/pep/tf_map.tsv
```

Command:

```bash
regutrace map-dbd-motifs \
  --tf-genes hf1-tf_genes.tsv \
  --proteins hf1-tfs.faa \
  --dbd-models config/dbd-models/plant-v1.json \
  --reference-map reference-map.tsv \
  --donor-dbd donor-dbd.tsv \
  --donor-dbd-faa donor-dbd.faa \
  --threshold-config config/dbd-thresholds/exploratory-v1.json \
  --allow-uncalibrated true \
  --reference-pack-version plant-jaspar-2024-rbh \
  --blastp <blastp> \
  --out hf1-dbd-map.tsv \
  --diagnostics-out hf1-dbd-diagnostics.tsv
```

The threshold config is explicitly `uncalibrated`; the command refuses to run
without `--allow-uncalibrated true`. Output rows are labelled
`mapping_confidence=exploratory_uncalibrated` and
`mapping_method=dbd_transfer_exploratory`. These are candidates, not validated
targets.

Result:

| item | value |
|---|---|
| catalogue TFs | 1,981 |
| TFs with any DBD candidate | 1,113 |
| candidate rows | 28,465 |
| unique donors among candidates | 565 |
| unique motifs among candidates | 565 |
| mapped with a single donor+motif | 59 |
| mapped with ambiguity | 1,054 |
| `NO_DBD_SEQUENCE` | 341 |
| `BELOW_DBD_THRESHOLD` | 203 |
| `NO_BLAST_HIT` | 114 |
| `NO_DONOR_DBD` | 89 |
| `NOT_IDENTIFIED` | 74 |
| `FAMILY_MISMATCH` | 47 |

The large ambiguity count is expected: the exploratory threshold admits many
same-family donors and the output preserves all of them instead of picking the
first one.

## Comparison with the strict HF1 map

The strict map has 260 mapped TFs. DBD candidates recover the exact
`(tf_id, reference_tf_id, motif_id)` strict triple for 237/260 = **91.2%**.

Per-family recovery of strict triples:

| family | strict rows | recovered | recovery |
|---|---:|---:|---:|
| ARF, B3, BES1, CAMTA, Dof, E2F/DP, EIL, FAR1, GATA, GRF, HD-ZIP, HSF, LBD, LFY, MIKC, MYB-related, NAC, Nin-like, RAV, SBP, Trihelix, VOZ, WRKY, bHLH, bZIP | 1–25 each | all | 100% |
| MYB | 31 | 30 | 96.8% |
| ERF | 27 | 25 | 92.6% |
| AP2 | 4 | 3 | 75.0% |
| BBR-BPC | 2 | 1 | 50.0% |
| TCP | 7 | 3 | 42.9% |
| C2H2 | 8 | 0 | 0% (abstained) |
| CPP | 2 | 0 | 0% (abstained) |
| ZF-HD | 4 | 0 | 0% (abstained) |

The highest-bitscore DBD donor equals the strict donor for 180/260 strict TFs
(69.2%). The remaining strict donors are still present among the candidates for
most TFs but are not ranked first, which is exactly why Phase 5 calibration and
benchmarking are required before this mode can be called balanced or strict.

Key HF1 site:

```text
TF:      Fta_6G0039340 (MYB)
strict:  Q9LTC4 -> MA0574.1
DBD:     best donor Q9LTC4, identity 82.83%, query coverage 100%,
         donor coverage 100%, E-value 3.12e-64
         (50 eligible same-family donors; ambiguity preserved)
```

This is a consistency check, not independent validation: the strict map is a
full-protein RBH result, and both methods share the same donor proteins and
motif library.

## Output contract

`map-dbd-motifs` writes the ten legacy `tf_map.tsv` columns plus:

```text
mapping_confidence, dbd_identity, dbd_bitscore, dbd_evalue,
query_dbd_coverage, donor_dbd_coverage, alignment_length,
query_dbd_start, query_dbd_end, donor_dbd_start, donor_dbd_end,
query_full_dbd_start, query_full_dbd_end,
donor_full_dbd_start, donor_full_dbd_end,
threshold_config_version, threshold_source,
architecture_status, ambiguity_status, reason_code,
reference_pack_version
```

The first ten columns remain readable by `io::read_tf_map` and the existing
scan/trace path. `map-dbd-motifs` is standalone in Phase 2; the pipeline still
uses the strict mapper by default.

## Tests

```bash
export PATH="$LOCALAPPDATA/w64devkit/w64devkit/bin:$PATH"
export CC="$LOCALAPPDATA/w64devkit/w64devkit/bin/gcc.exe"
export CXX="$LOCALAPPDATA/w64devkit/w64devkit/bin/g++.exe"
export RUSTFLAGS="-C link-self-contained=yes"
cargo +stable-x86_64-pc-windows-gnu test --release --locked
```

Result:

```text
38 CLI unit tests passed (Phase 2 snapshot; current 40)
6 CLI integration tests passed (current 7)
69 core unit tests passed
0 failed
```

Phase 2-specific tests:

- `dbd_candidates_apply_family_and_threshold_filters` — family gate, below
  threshold, family mismatch, no-hit reasons;
- `dbd_candidates_keep_best_alignment_per_donor` — one candidate per donor;
- `phase2_dbd_extract_and_mapping_are_auditable` — `dbd-extract` and
  `map-dbd-motifs` end-to-end with a synthetic fixture and precomputed
  alignments (no external BLAST dependency in tests).

No strict mapping or scan regression was introduced.

## Scientific limitations

1. **Thresholds are not calibrated.** `plant-dbd-exploratory-v1` is an
   engineering placeholder (50% identity, 70% coverage, 1e-5). It must not be
   presented as a validated DBD transfer threshold. The tool requires
   `--allow-uncalibrated true` and labels every row.
2. **DBD boundaries are curated, not learned.** The config abstains where the
   domain architecture is variable or composite. C2H2/C3H arrays are the main
   coverage loss.
3. **One DBD representative per gene.** Isoforms are reduced to the
   best-scoring architecture; isoform-level DBD differences are not reported.
4. **RAV spans two domains.** The extracted region runs from the first AP2
   residue to the last B3 residue, including the linker. This is intentionally
   labelled as `full_architecture`, but it is not a minimal DBD.
5. **No independent benchmark.** The 91.2% strict-map recovery is a consistency
   check. It does not establish that DBD transfer is more accurate than strict
   RBH or that any candidate is a real TF–target relationship.
6. **External sources.** CIS-BP and PlantTFDB are accessible for local
   analysis, but redistribution remains unverified; no external matrices were
   added to this phase.

## Compatibility and migration

- No changes to `tf_map.tsv` semantics, `evidence.tsv`, `evidence.db`, or the
  pipeline.
- New files and configs are additive.
- A future pack should ship `donor-dbd.tsv` / `donor-dbd.faa` and the DBD model
  version; until then, they are generated locally from the pack proteins.
- The output has extra columns; strict consumers that enforce an exact column
  count must be updated, but header-based readers are unaffected.

## Next steps

1. Add C2H2/C3H array extraction with an explicit array boundary rule and a
   separate model version.
2. Calibrate DBD identity/coverage on a donor-held-out benchmark; do not reuse
   full-protein RBH thresholds.
3. Add a balanced mode that ranks candidates and requires a minimum margin or
   motif agreement; keep `exploratory` available.
4. Integrate the DBD candidate table into trace/dashboard only after Phase 4
   evidence-tier UX is defined.
5. Build JASPAR 2026 and PlantTFDB adapters only after provenance and
   redistribution terms are recorded.
