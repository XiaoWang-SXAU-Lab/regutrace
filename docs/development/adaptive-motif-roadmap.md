# Adaptive motif mapping and evidence-aware target prediction

Status: **Phase 0 complete; Phase 1 legacy registry and JASPAR 2026 candidate
pack complete, external-source integration partly blocked; Phase 2 exploratory
DBD transfer complete; Phase 2.5 calibration complete with two calibrated
families; Phase 2.6 mapping-strategy ablation complete; Phase 3 scan
architecture, statistics and five HF1 benchmarks complete** (2026-11-25). This
document records the real architecture, audit findings, planned changes and
acceptance criteria for the ReguTrace upgrade. It is intentionally updated after
every phase.

## Scope

Goal: improve motif coverage for non-model plants while preserving the
distinction between:

1. TF identification;
2. TF–motif mapping;
3. motif–DNA occurrence;
4. TF–target candidate relationships.

The current strict RBH workflow remains the reference baseline. New methods must
be additive, versioned and independently labelled. A mapped motif never upgrades
a `level D` promoter hit to experimental evidence.

## Actual call chain

| Stage | Code | Output |
|---|---|---|
| CLI dispatch | `crates/regutrace-cli/src/main.rs` (`pipeline`, `map-tf-motifs`, `identify-tfs`, `scan`, `trace`) | CLI |
| TF identification | `crates/regutrace-cli/src/plant_tfs.rs` (`run`, `prepare_fasta`, `write_catalogue`, `parse_gff_gene_map`) | `tf_genes.tsv` |
| Family/domain classification | `crates/regutrace-core/src/tf_classification.rs`, `hmmer.rs` | classification + domain evidence |
| TF–motif mapping | `crates/regutrace-cli/src/motif_transfer.rs` (`automatic_hits`, `best_hits_by`, `write_local_domains`) | `tf_map.tsv` |
| Pipeline orchestration | `crates/regutrace-cli/src/pipeline.rs` (`run`, `output_paths`, report assembly) | `pipeline-report.json` |
| Promoter extraction | `crates/regutrace-core/src/promoter.rs`, `gff.rs` | `promoters.fa` |
| PWM scan | `crates/regutrace-core/src/pwm.rs`, `scan.rs` | `evidence.tsv` |
| Trace views | `crates/regutrace-cli/src/trace.rs` | `tf_motif_trace.tsv`, `motif_ambiguity.tsv`, `target_candidate_ties.tsv` |
| Persistence | `crates/regutrace-core/src/store.rs` | `evidence.db`, `prior` |
| Dashboard | `wgcna-dashboard/backend/internal/regutrace`, `frontend/src/views/ToolkitRegutraceView.vue` | job API/UI |
| Plant reference pack | `plant-jaspar-2024-rbh` (`manifest.json`, `reference-proteins.faa`, `reference-map.tsv`) | 615 donors, one representative JASPAR matrix each |

## Current strict algorithm (legacy baseline)

1. `plant_tfs` identifies TF genes from proteins/CDS and writes a complete
   `tf_genes.tsv`; ambiguous genes are retained as ambiguous.
2. For the new plant pack, `pipeline` passes `--reference-proteins` to
   `map-tf-motifs`. `automatic_hits` writes all protein records for **identified**
   genes, runs BLASTP local→donor and donor→local, and maps only if:
   - the local gene is `identified`;
   - local and donor family labels match (semicolon-separated labels accepted);
   - the local→donor best hit contains the donor;
   - the donor→local reverse best hit contains the local gene;
   - the best alignment passes identity ≥35%, query/donor coverage ≥60%, and
     E-value ≤1e-5.
3. Ties are retained as multiple rows and labelled `*_tied_candidate`.
4. `scan` filters single-domain fallback mappings unless explicitly enabled.
5. Promoter PWM hits are `evidence_level=D`; trace tables preserve donor/motif
   provenance.

The exact thresholds are implementation defaults in
`motif_transfer.rs::best_hits_by`; they are **not** yet a versioned,
family-specific scientific calibration.

## Phase 0 audit findings

### P0-A1 — best-hit threshold is applied only to the top-scoring alignment

`best_hits_by` keeps only the highest BLAST bitscore per query and then records
whether that best row passes identity/coverage/E-value. If the top raw hit is
ineligible, a lower-scoring eligible hit is not considered. This can hide a
valid mapping and makes the unmapped reason ambiguous. The legacy behaviour is
kept for now; diagnostics must expose it rather than silently changing results.

### P0-A2 — family compatibility is a hard gate, not a scored candidate property

Whole-protein BLAST rows are filtered by family before ranking. A local TF with
no family-compatible donor is indistinguishable from a TF with no BLAST hit in
the current report. Reason codes must distinguish:

- no family donor/motif in the reference map;
- BLAST hits but none in a compatible family;
- best raw hit below thresholds;
- best raw hit not reciprocal.

### P0-A3 — ambiguous genes are excluded from the whole-protein query

`automatic_hits` inserts only `gene.is_identified()` into `local_index` for the
`--reference-proteins` branch. Ambiguous TF genes are therefore not queried,
while the legacy DBD path permits ambiguous genes to compete but rejects them at
the final mapping step. This is an intentional-looking inconsistency and must be
documented and tested.

### P0-A4 — full-length protein RBH is not DBD-specific transfer

The current plant workflow uses the full protein for BLASTP. It does not use the
HMMER DBD coordinates, DBD type or architecture for ranking. The strict result
must remain available as `Legacy strict`; DBD-based transfer belongs in Phase 2.

### P0-A5 — multi-isoform and domain extraction are lossy by design

`write_local_domains` chooses one domain (or one combined architecture) per gene
by domain score; the whole-protein path keeps each protein sequence and lets BLAST
pick the best. The catalogue can contain multiple isoforms but mapping does not
produce an isoform-level audit. Phase 0 should report the gene-level primary
reason without inventing isoform-level claims.

### P0-A6 — no per-TF / per-family unmapped diagnostics

`pipeline-report.json` currently reports only:
- `tf_catalogue` counts;
- `tf_map.tf_without_motif.count` and 100 examples.
It cannot answer why each of the 1,648 HF1 TFs is unmapped. There is no
per-family donor/motif count or reason breakdown.

### P0-A7 — reference map has one row per donor but no registry provenance

The pack manifest records checksums and a version, but not per-donor source,
experimental method, DBD coordinates, motif provenance or license. Multiple
motifs per donor cannot be represented by the current `reference-proteins`
contract (it rejects duplicate donor IDs). This is Phase 1 work.

### P0-A8 — family-label conflicts are not preserved

The HF1 audit case `Fta_3G0045350` is classified `MYB` by the domain rules while
its HF1/PlantTFDB annotation is `MYB-related`/`DIVARICATA-like` and is flagged
`different_label_review`. Mapping treats the catalogue family as authoritative;
the conflict is not represented in `tf_map` diagnostics. Phase 0 must expose it
as an audit reason, not silently resolve it.

### P0-A9 — scanning statistics are not calibrated

The current scan uses a uniform background, 1 kb upstream, per-motif p≤1e-4, no
genome-wide FDR and a target union. These are candidate filters; Phase 0 must
record them as such and not present 37,407/37,415 as a biological validation.

### P0-A10 — frontend currently derives family mapped counts from summary data

The dashboard has an existing family chart (`tfFamilyRows`) and target/TF detail
views, but does not consume a per-TF mapping reason or a per-family diagnostic
contract. Phase 4 can render Phase 0 diagnostics without changing the visual
language.

## Phase 0 implementation

### Outputs

Add two optional, backward-compatible mapping diagnostics:

- `mapping-diagnostics.tsv` — one row per catalogue TF:
  `tf_id`, `tf_symbol`, `family`, `status`, `mapping_status`,
  `primary_reason`, `all_reasons`, `mapping_method`, `motif_id`,
  `reference_tf_id`, `forward_best_donor`, `forward_best_score`,
  `forward_best_eligible`, `reverse_best_recipients`, `family_donor_count`,
  `family_motif_count`, `notes`.
- `mapping-diagnostics.family.tsv` — one row per exact catalogue family:
  catalogue/identified/ambiguous/rejected counts, mapped TFs and rows,
  unmapped identified, primary-reason counts, unique donors and motifs.

Reason vocabulary (actual logic, not a score):

| code | meaning |
|---|---|
| `MAPPED` | one or more strict reciprocal mappings |
| `MAPPED_TIED` | strict mapping but tied donor/recipient |
| `NOT_IDENTIFIED` | ambiguous/rejected catalogue row not eligible for final mapping |
| `NO_DBD_SEQUENCE` | DBD path could not extract a configured DBD sequence |
| `NO_REFERENCE_DONOR` | no reference-map donor in the TF's family |
| `NO_BLAST_HIT` | no local→donor alignment survived the family filter |
| `BELOW_IDENTITY_OR_COVERAGE` | best family-compatible alignment failed thresholds |
| `NOT_RECIPROCAL_BEST` | forward best donor does not return this gene as reverse best |
| `FAMILY_MISMATCH` | mapping candidates exist but family labels do not overlap |
| `MULTIPLE_BEST_DONORS` | multiple equally ranked donors/recipients |
| `UNRESOLVED` | no rule explains the state; must be investigated |

`all_reasons` can contain several codes; `primary_reason` is the highest-priority
actual blocker. The summary counts only the primary reason, so rows are not
double-counted; all reasons remain auditable in the per-TF table.

### Pipeline integration

- `map-tf-motifs` accepts optional `--diagnostics-out` and
  `--diagnostics-family-out`; existing invocations and outputs are unchanged.
- `pipeline` writes both files and adds a `mapping_diagnostics` block to
  `pipeline-report.json` with paths and per-family rows.
- The strict `tf_map.tsv` is unchanged; no reason code changes mapping.

### Tests

- A synthetic explicit forward/reverse fixture exercises:
  mapped, mapped-tied, no donor, no BLAST hit, below threshold, non-reciprocal,
  family mismatch, ambiguous excluded, and no-DBD path.
- Assert per-family counts sum to the catalogue total without double counting.
- Assert diagnostics are deterministic and do not alter `tf_map.tsv`.
- Run `cargo test --release --locked` and a small real-pack smoke test when
  available.
- HF1 is a **reference case**, not an automated fixture: do not claim a full
  reproduction unless the full input and pack are present.

## Phase 1 — reference registry and pack builder

Implemented in the legacy scope (see [phase1-report.md](phase1-report.md)):

- versioned `RegistryRecord` schema with source, method, donor protein, species,
  taxon ID, DBD coordinates, license, checksum and QC status;
- `registry inspect` for legacy packs and registry TSVs;
- per-donor species/tax ID from the legacy reference map;
- per-donor `sha256:` donor-sequence checksums with mismatch validation;
- manifest checksum and file-presence verification;
- duplicate-ID, missing-sequence and DBD-coordinate validation;
- motif-library check (`--motif-library`): registry motifs absent from the
  library fail; unused library motifs and exact normalized-matrix duplicate
  groups are reported;
- per-source, per-family and per-species summary JSON;
- end-to-end tests with a synthetic legacy pack.

Audited sources:

- JASPAR 2024 pack library is the official JASPAR 2024 CORE plants file
  (SHA-256 exact match, 805/805 IDs);
- JASPAR 2026 release 11 (2025-12-31) is CC BY 4.0 and has 927 plant CORE
  matrices. A candidate `plant-jaspar-2026-rbh` pack was built (713 donors,
  713 representative motifs, 34 families) and stored in the sync analysis
  folder with a NOTICE. HF1 strict mapping: 271 vs 260 mapped TFs, 250 common,
  21 new, 10 lost, 23 changed. The 2024 pack remains the default until
  target-level benchmarking; no external motifs were committed to the code
  repository;
- CIS-BP 3.10 archives are accessible by URL but the 240 MB GET was throttled to
  timeout from this network; the RSAT mirror works, but its LICENSE.txt does not
  grant redistribution. CIS-BP remains redistribution-blocked;
- PlantTFDB experimental data is accessible (674 MEME motifs, 675-row info
  table) but no explicit data license was found; redistribution unverified;
- PlantRegMap, Plant Cistrome and PlantPAN are not yet audited.

Not implemented / blocked:

- rebuilding the pack from the registry (the builder script depends on
  `donor-audit.tsv` and `jaspar-tfs.faa`, absent from the repository);
- external-source adapters;
- per-donor DBD coordinates: `reference-domains.faa` uses Arabidopsis AGI IDs
  and has zero overlap with the 615 donor accessions. Phase 2 must extract donor
  DBDs with HMMER.

## Phase 2 — adaptive mapping

Implemented (see [phase2-report.md](phase2-report.md)):

- Tier 1 legacy strict RBH remains the default and is unchanged.
- `dbd-extract` extracts one DBD per gene from a catalogue using a versioned
  family-to-model config; variable/composite families abstain.
- `map-dbd-motifs` aligns query DBDs to donor DBDs, applies a family gate and
  versioned threshold config, retains all eligible donors, and writes explicit
  provenance fields:
  `mapping_confidence`, `dbd_identity`, `dbd_bitscore`, `dbd_evalue`,
  `query_dbd_coverage`, `donor_dbd_coverage`, `threshold_config_version`,
  `threshold_source`, `architecture_status`, `ambiguity_status`, `reason_code`,
  `reference_pack_version`.
- `config/dbd-models/plant-v1.json`: 45 supported families, 10 abstained.
- `config/dbd-thresholds/exploratory-v1.json`: uncalibrated; requires
  `--allow-uncalibrated true`.
- Real pack: 585/615 donor DBDs extracted; 30 abstained (C2H2, CPP, ZF-HD).
- HF1 exploratory run: 1,113 TFs with candidates, 28,465 rows, 91.2% recovery of
  the strict donor+motif triples; 180/260 strict donors rank first by bitscore.
- All outputs are labelled exploratory; nothing is mixed into `tf_map.tsv`.

Still required before a balanced/strict DBD mode:

- family-specific DBD identity/coverage calibration with donor-held-out
  benchmark;
- C2H2/C3H array extraction and a separate model version;
- candidate ranking/margin or motif-agreement rule to control ambiguity;
- pipeline/dashboard integration after Phase 4 evidence-tier UX;
- independent biological benchmark (DAP-seq/ChIP-seq) in Phase 5.

Modes: `legacy_strict` remains default; `dbd_exploratory` is available as a
standalone command. `balanced` and `strict` DBD modes are not implemented.

## Phase 2.5 — DBD transfer calibration

Implemented between Phase 2 and Phase 3 (see
[phase2.5-report.md](phase2.5-report.md),
[dbd-calibration-methodology.md](dbd-calibration-methodology.md) and
[phase2.5-decisions.md](phase2.5-decisions.md)):

- 596-donor / 619-motif calibration dataset; 2024 and 2026 releases merged by
  donor and motif; exact PWM groups; experimental-method filter; explicit
  exclusions.
- TOMTOM-style PWM correlation (offset and reverse-complement aware, ≥50%
  overlap), with 0.80 primary and 0.90 sensitivity labels.
- 13,848 within-family donor pairs; DBD identity predicts motif similarity
  better than full-length protein identity overall (AUC 0.696 vs 0.608), but not
  in every family.
- 174 specificity groups at the 0.90 label (619 motifs; 61 singletons).
- Nested species-grouped cross-validation; only `HD-ZIP` (75% identity / 0.90
  coverage) and `TCP` (70% / 0.50) met the promotion criteria.
- 30 families abstain; no universal threshold and no coverage target.
- Calibrated mapping statuses in `map-dbd-motifs`:
  `dbd_calibrated`, `dbd_ambiguous`, `dbd_exploratory`,
  `no_supported_motif`; reciprocal DBD-best and specificity-group ambiguity are
  required for a calibrated call.
- HF1: calibrated DBD maps 19 TFs / 57 rows; strict 2026 + calibrated maps 280
  TFs, adding 9 HD-ZIP TFs. Ambiguous DBD TFs are excluded from the combined set.
- 29 short motifs have an attainable p-value floor above `1e-4`; the cutoff was
  not silently relaxed.

Not implemented in Phase 2.5:

- C2H2/C3H array extraction;
- DNA-contact-residue features;
- independent DAP-seq/ChIP-seq benchmark;
- Phase 3 scanning of calibrated associations.

## Phase 2.6 — mapping-strategy ablation

Implemented after Phase 2.5 (see
[phase2.6-report.md](phase2.6-report.md),
[phase2.6-cisbp-method-review.md](phase2.6-cisbp-method-review.md) and
[phase2.6-decisions.md](phase2.6-decisions.md)):

- Compared A full-length strict RBH, B DBD reciprocal, C DBD non-reciprocal
  global threshold, D DBD non-reciprocal family thresholds, and E
  specificity-group-aware inference on the same 596-donor calibration set and
  outer folds.
- Species-held-out precision/coverage: A 0.878/0.233, B 0.933/0.151,
  C 0.767/0.072, D 0.846/0.272, E 0.864/0.074.
- Reciprocal DBD (B) raises precision but lowers coverage; non-reciprocal DBD
  (D) recovers more correct calls with lower precision and does not meet the
  0.90 target under homology-cluster-held-out evaluation.
- Validation setting changes the estimate: donor-held-out is optimistic;
  species- and homology-cluster-held-out are lower. The three settings are not
  equivalent.
- Pairwise AUC CIs are pseudoreplicated at pair level: 53% of evaluable pairs
  are within a 70% DBD homology cluster; cluster bootstrap CI width is
  ~5–6× the pair bootstrap width.
- No Phase 2.5 calibrated config was overwritten; no non-reciprocal method was
  promoted; no UI or genome-wide scan was performed.

## Phase 3 — scan/storage/statistics

Reports:
[phase3-report.md](phase3-report.md),
[phase3-scanner-audit.md](phase3-scanner-audit.md),
[phase3-statistical-methodology.md](phase3-statistical-methodology.md),
[phase3-hf1-comparison.md](phase3-hf1-comparison.md).

- implemented exact normalized-matrix hashing, opt-in exact PWM reuse and a
  persistent content-keyed scan cache;
- implemented a unique site table, promoter summary and testability table;
  shared motif sites are counted once and TF evidence rows remain compatible;
- implemented a strand-symmetric empirical zero-order promoter background with
  provenance; uniform remains the default;
- exact scanner p-values and attainable floors are reported; no site/promoter
  q-values are claimed; exploratory Fisher/BH enrichment is available;
- legacy and optimized HF1 evidence are full-row identical; peak memory fell
  from 3,308 to 454 MiB; warm-cache scans are 4–5× faster;
- FIMO comparison and normalized SQLite scan/site/aggregation tables remain
  deferred.

## Phase 4 — UI/evidence UX

- overview: identified / eligible / mapped / scanned / with-hit / unmapped reasons;
- per-TF diagnostics and donor candidates;
- TF-specific vs family exploratory filtering;
- empty-state explanations: `No reliable motif assigned` vs
  `Scanned; no hits at current threshold`;
- prior export only from strict/balanced tiers.

## Phase 5 — benchmarks and ablation

- leave-one-TF-out, homology-group and species splits; remove database
  duplicates and threshold leakage;
- DAP-seq/ChIP-seq target benchmarks with explicit label limitations;
- ablation A–E and HF1 re-run only when full inputs are available;
- no fabricated metrics.

## Phase 6–7 — integration/release

Only after Phase 0–5: motif enrichment + expression evidence in rsWGCNA,
versioned export metadata, packaging, migration, documentation, licensing and
cross-platform release checks.

## Compatibility and risks

- `tf_map.tsv`, `evidence.tsv`, `tf_genes.tsv` and `pipeline-report.json` must
  remain readable by the dashboard and prior exporter.
- New diagnostics are additive; old jobs without them remain valid.
- Strict mapping semantics must not change in Phase 0.
- Full-protein RBH and DBD transfer must not be conflated in labels.
- HF1 exact counts are environment-dependent (macOS vs Windows toolchain); do not
  treat small evidence deltas as logic regressions without an input manifest.
- Reference pack redistribution/licensing is a blocker for Phase 1 sources until
  audited.

## Acceptance status

- [x] Repository and actual call chain audited.
- [x] Baseline constraints and known failure modes recorded.
- [x] Per-TF reason diagnostics implemented.
- [x] Per-family summary implemented.
- [x] Synthetic fixture and tests pass (Phase 0 baseline: 31 CLI + 4 integration
  + 69 core; current suite: 40 CLI + 7 integration + 69 core).
- [x] Strict output backward compatibility checked on the synthetic fixture.
- [x] Phase 0 report updated with actual commands/results
  ([phase0-report.md](phase0-report.md)).
- [x] Phase 1 legacy registry implemented and validated on the 615-donor pack.
- [x] Phase 1 source/license audit recorded
  ([phase1-report.md](phase1-report.md)).
- [x] Phase 1 pack rebuild demonstrated with the local donor audit inputs.
- [x] Phase 1 JASPAR 2026 candidate pack built and validated.
- [ ] Phase 1 CIS-BP / PlantTFDB redistribution-compatible adapters (blocked).
- [x] Phase 2 DBD extraction and exploratory transfer implemented.
- [x] Phase 2 real-pack donor DBD extraction validated (585/615; 0 mismatches).
- [x] Phase 2 HF1 consistency check recorded (91.2% strict triple recovery).
- [x] Phase 2.5 calibration dataset, pairwise analysis and specificity groups
  implemented ([phase2.5-report.md](phase2.5-report.md)).
- [x] Phase 2.5 family-specific thresholds calibrated for HD-ZIP and TCP; 30
  families explicitly abstain.
- [x] Phase 2.5 HF1 reassessment and short-motif p-value audit
  ([phase2.5-hf1-comparison.md](phase2.5-hf1-comparison.md)).
- [x] Phase 2.6 mapping-strategy ablation and validation-setting audit
  ([phase2.6-report.md](phase2.6-report.md)).
- [x] Phase 3 scanner audit, exact PWM reuse, site/testability layers,
  empirical background, synthetic calibration and five HF1 runs
  ([phase3-report.md](phase3-report.md)).
- [x] Phase 3 legacy/optimized evidence equivalence and memory reduction.
- [x] Phase 3 exploratory motif enrichment interface and tests.
- [ ] Phase 3 normalized SQLite scan/site/aggregation schema (deferred).
- [ ] Phase 3 site/promoter q-values with a dependence-aware null (deferred).
- [ ] Phase 3 FIMO comparison (FIMO not installed).
- [ ] Phase 2 calibrated balanced/strict genome-wide DBD mode (requires Phase 5
  benchmark and Phase 3 scan integration; Phase 3 scan integration is available,
  Phase 5 benchmark is not).
