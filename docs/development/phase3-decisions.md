# Phase 3 decisions

Date: 2026-11-25

## D1. The legacy scan path is the default

`scan` uses the original motif-ID loop, in-memory sorting and exact legacy
evidence rows unless a Phase 3 option is supplied. This preserves dashboard,
`pairs` and `prior` behavior. The Phase 3 path is opt-in via `--reuse-exact-pwm`,
`--scan-cache`, `--sites-out`, `--motif-testability-out`,
`--promoter-summary-out`, `--scan-manifest` or `--background-mode empirical`.

## D2. Exact reuse means bit-identical normalized matrices

`Pwm::digest()` hashes the normalized probability bits, not the accession, name,
family or a similarity cluster. There is no similarity threshold in the reuse
layer. Two matrices that differ in one floating-point bit are scanned separately.

This is deliberate: merging similar matrices is a mapping decision, not a scan
optimization. The Phase 2.5/2.6 similarity work remains at the mapping layer.

## D3. Cache keys include every scan-determining input

The key covers algorithm version, score scale, background values, p-value,
min-score, strand mode, promoter-set hash and PWM digest. A changed p-value,
background or promoter set invalidates the cache. The p-value stays in the key
because the cache stores hits, not the full score distribution.

No similarity, family or accession is in the key. The cache is content-addressed,
not annotation-addressed.

## D4. Site evidence and TF evidence are separate

The site table has one row per unique PWM group and site. The evidence table has
one row per TF, motif ID and site to preserve legacy behavior. The two counts are
reported separately:

- `unique_site_hits`;
- `motif_id_hits`;
- `evidence_rows`.

A shared motif site is not independent evidence for every TF. The dashboard and
downstream statistics should use the site table for site counts and the evidence
table for TF–target associations.

## D5. Empirical background is a scientific option, not a default

`--background-mode empirical` uses a strand-symmetric zero-order model with a
pseudocount of 1 per base, estimated from the supplied promoters. It is not the
default because:

- it changes the null and the hit set materially (16.9% fewer sites on HF1);
- it is estimated from the same promoters that are scanned, so p-values are
  conditional on the estimate;
- it is still zero-order and does not model dinucleotides or GC windows.

Uniform remains the default for comparability with all prior phases. The
manifest records the model and raw counts so a run can be audited.

## D6. No site or promoter q-values in Phase 3

The scanner p-value is exact for one pre-specified window under the model.
Overlapping windows, two strands, multiple motifs and shared PWMs induce
dependence. BH over all sites would not control the FDR for a biological
question, and BH over only significant sites would be invalid.

Phase 3 therefore reports:

- exact site p-values;
- attainable p-value floors;
- descriptive promoter summaries;
- exploratory Fisher/BH enrichment.

It does not report a genome-wide site FDR or a promoter q-value. The required
future work is a sequence-based null preserving length/GC/dinucleotide
composition and a block resampling or dependence-aware procedure.

## D7. Testability is explicit and non-destructive

A motif whose exact floor is above the requested p-value is marked
`threshold_unattainable` and skipped. The requested threshold is never lowered
silently. This affects 28/257 (10.9%) of 2024 mapped motifs and 29/268 (10.8%) of
2026 mapped motifs under uniform; 19/268 (7.1%) under the empirical background.

## D8. Enrichment is exploratory

`regutrace enrich` implements Fisher/BH on a binary promoter-hit matrix over the
scanned universe. It is tested and documented, but not calibrated for promoter
length, GC, motif similarity or gene-set overlap. It is not used to promote any
mapping or evidence tier.

## D9. The SQLite compatibility projection is retained

`evidence.db` keeps the `hits`/`pairs` schema. Normalized scan/site/aggregation
tables are not added yet. The new layers are TSV/JSON outputs. This avoids
breaking `prior`, `validate-db` and the dashboard while the Phase 3 interfaces
settle.

## D10. Phase 3 does not change mapping or evidence tiers

No Phase 2.5 threshold was changed. No non-reciprocal DBD method was promoted.
No `dbd_exploratory` row was added to the HF1 run. The 9 calibrated HD-ZIP TFs in
Run 5 are `dbd_calibrated` and remain evidence level D.

## D11. HF1 Run 5 uses a provenance-preserving map

The Phase 2.5 `hf1-strict-plus-calibrated.tsv` summary lacks
`mapping_method`/`domain_mode`, so using it directly changes evidence `notes` for
strict rows. Phase 3 builds
`hf1-strict-plus-calibrated-with-provenance.tsv` by taking the 271 strict rows
from the Phase 3 map and adding the 9 calibrated HD-ZIP rows from
`hf1-dbd-calibrated.tsv`. This makes the strict subset exactly equal to Run 3 and
keeps calibrated provenance in the notes.

## D12. No Phase 2.7 or threshold tuning

Phase 3 stops at scan architecture, statistics and HF1 benchmarking. It does not
start Phase 2.7 or continue tuning family thresholds on the Phase 2.5 calibration
set. That work requires the Phase 5 external benchmark.

## D13. No commit, push, merge or installer rebuild

The working tree is uncommitted. No installer was rebuilt. No branch was
switched or merged. The sync workspace contains data and documentation only.
