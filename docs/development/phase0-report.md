# Phase 0 report — audit, baseline and mapping diagnostics

Date: 2026-11-24  
Repository: `regutrace`  
Base HEAD: `9ac4837` (plus this working tree; not committed)

## Scope delivered

Phase 0 is an additive audit/diagnostic layer. The strict mapping output
(`tf_map.tsv`) is not changed.

Implemented:

- per-TF `mapping-diagnostics.tsv`;
- per-family `mapping-diagnostics.family.tsv`;
- reason vocabulary documented in
  [adaptive-motif-roadmap.md](adaptive-motif-roadmap.md);
- `map-tf-motifs` optional flags:
  - `--diagnostics-out <tsv>`
  - `--diagnostics-family-out <tsv>`
- `pipeline` writes both files and exposes them in
  `pipeline-report.json.mapping_diagnostics`;
- a synthetic fixture covering mapped, tied, no donor, no BLAST, below
  threshold, non-reciprocal, family mismatch, ambiguous and rejected cases.

## Actual code path audited

- `crates/regutrace-cli/src/plant_tfs.rs` — TF identification and catalogue.
- `crates/regutrace-cli/src/motif_transfer.rs` — family filter, RBH, tie logic,
  local-domain extraction.
- `crates/regutrace-cli/src/pipeline.rs` — orchestration and report assembly.
- `crates/regutrace-core/src/io.rs` — `tf_genes`/`tf_map` schema.
- `crates/regutrace-core/src/scan.rs`, `pwm.rs`, `promoter.rs` — scanning.
- `crates/regutrace-cli/src/trace.rs` — trace outputs.
- Dashboard: `wgcna-dashboard/backend/internal/regutrace` and
  `frontend/src/views/ToolkitRegutraceView.vue` (existing family/TF views; no
  Phase 0 frontend change).

The current plant pack path uses `--reference-proteins` and full-protein
reciprocal BLASTP; the legacy DBD path remains available for old packs.

## Important audit findings

1. The best hit is selected by raw BLAST bitscore first; identity/coverage/
   E-value are applied to that best row. A lower eligible hit is not substituted.
   This is preserved for strict compatibility and should be challenged only behind
   a new mode.
2. Whole-protein family filtering happens before ranking; the new diagnostics
   distinguish no alignment from no family-compatible alignment.
3. Ambiguous genes are not inserted into the whole-protein query; the legacy DBD
   path behaves differently. The diagnostic output marks both as
   `AMBIGUOUS_EXCLUDED`.
4. No per-TF reason was previously persisted; only a count and 100 examples of
   `tf_without_motif` existed.
5. The 615-donor pack has one representative matrix per donor; `reference-proteins`
   rejects duplicate donor IDs. Multi-motif donor provenance is Phase 1.
6. Family label conflicts (for example `Fta_3G0045350`: classifier `MYB`,
   PlantTFDB `MYB-related`/DIVARICATA-like) are not resolved by the mapper and
   must remain visible in diagnostics.

## Diagnostics contract

Per-TF fields:

```text
tf_id, tf_symbol, family, status, mapping_status, primary_reason, all_reasons,
mapping_method, motif_id, reference_tf_id, forward_best_donor,
forward_best_score, forward_best_eligible, reverse_best_recipients,
family_donor_count, family_motif_count, notes
```

`primary_reason` is the highest-priority actual blocker; `all_reasons` retains
all observed reasons. Family summary counts only primary reasons so rows are not
double-counted.

Reason codes currently emitted:

```text
MAPPED
MAPPED_TIED
AMBIGUOUS_EXCLUDED
REJECTED_EXCLUDED
NO_DBD_SEQUENCE
NO_REFERENCE_DONOR
NO_BLAST_HIT
FAMILY_MISMATCH
BELOW_IDENTITY_OR_COVERAGE
NOT_RECIPROCAL_BEST
MULTIPLE_BEST_DONORS
UNRESOLVED
```

## Synthetic fixture

Location:

```text
crates/regutrace-cli/tests/fixtures/phase0/
  tf_genes.tsv
  reference-map.tsv
  forward.tsv
  reverse.tsv
```

Expected primary reasons verified by
`phase0_mapping_diagnostics_are_auditable`:

| gene | primary reason |
|---|---|
| `G_MAPPED` | `MAPPED` |
| `G_TIED` | `MAPPED_TIED` |
| `G_NODONOR` | `NO_REFERENCE_DONOR` |
| `G_NOBLAST` | `NO_BLAST_HIT` |
| `G_BELOW` | `BELOW_IDENTITY_OR_COVERAGE` |
| `G_NONRECIP` | `NOT_RECIPROCAL_BEST` |
| `G_FAMILY` | `NO_REFERENCE_DONOR` with `FAMILY_MISMATCH` also present |
| `G_AMB` | `AMBIGUOUS_EXCLUDED` |
| `G_REJ` | `REJECTED_EXCLUDED` |

Family consistency check: family A has 7 catalogue rows, 5 identified, 2 mapped
TFs, 3 unmapped identified; primary-reason counts sum to 3.

## Tests run

Command:

```bash
export PATH="$LOCALAPPDATA/w64devkit/w64devkit/bin:$PATH"
export CC="$LOCALAPPDATA/w64devkit/w64devkit/bin/gcc.exe"
export CXX="$LOCALAPPDATA/w64devkit/w64devkit/bin/g++.exe"
export RUSTFLAGS="-C link-self-contained=yes"
cargo +stable-x86_64-pc-windows-gnu test --release --locked
```

Result:

```text
31 CLI unit tests passed
4 CLI integration tests passed
69 core unit tests passed
0 failed
```

The first attempt without the GNU toolchain failed in `libsqlite3-sys` because
`gcc.exe` and `-lgcc_eh` were unavailable; this is an environment issue, not a
code failure. The command above is the reproducible one on this machine.

## Pipeline smoke test

A small synthetic genome/GFF/protein was run through the real
`regutrace pipeline` with `plant-jaspar-2024-rbh` and bundled HMMER/BLASTP:

```text
exit 0
TF catalogue: 0 identified / 0 total
TF–motif: 0 rows
scan: 0 motif hits / 0 evidence rows
mapping-diagnostics.tsv: header only
mapping-diagnostics.family.tsv: header only
pipeline-report.json.mapping_diagnostics present
```

This verifies the pipeline writes the new artifacts without changing the empty
result path. It is not a biological validation and does not reproduce HF1.

## Compatibility

- Existing `tf_map.tsv` columns and row selection logic are unchanged.
- New diagnostics are optional for standalone `map-tf-motifs`.
- Pipeline output directory now includes two additional files; clients that
  enforce an exact output allowlist must be updated.
- Old jobs and packs without diagnostics remain readable.
- No database schema migration was needed.

## Remaining Phase 0 work / blockers

- Full HF1 strict regression is not rerun because the complete HF1 input and
  pack are not all in this repository; the report and the synthetic fixture are
  the available baseline.
- `motif missing from the pack library` is not yet a per-TF reason. The mapper
  does not receive the motif library; this belongs to Phase 1/3 integration.
- A best raw hit that is ineligible while a lower eligible hit exists is
  recorded as `BELOW_IDENTITY_OR_COVERAGE`; the lower eligible score is audited
  internally but not yet a separate reason/code.
- Multi-family catalogue rows are grouped by their exact semicolon-joined family
  string in the family summary; a later registry should define unambiguous
  family membership.

## Next phase

Phase 1: reference registry and pack builder. External JASPAR 2026/CIS-BP data,
licenses and donor accessions must be verified before download/redistribution.
Until then, registry schema and synthetic adapter tests can proceed, but no
source coverage or benchmark claim is permitted.
