# Phase 4 user workflows

**Date:** 2026-11-25

## 1. Product position

ReguTrace is a **standalone toolkit for genome-wide transcription factor identification and regulatory background construction**.

It does not require rsWGCNA, expression data, traits, or a WGCNA design. It builds a reusable regulatory background once for a genome/reference pack; rsWGCNA projects consume that background later.

Responsibility boundary:

| ReguTrace owns | rsWGCNA owns |
|---|---|
| genome-wide TF identification | expression-specific co-expression networks |
| TF–motif mapping and provenance | module construction |
| promoter motif scanning and testability | module–trait association |
| candidate TF–target prior (level D) | regulator prioritization within a project |
| reference-pack and scan provenance | project-specific filtering/ranking |

ReguTrace does not filter or rank by expression, module membership, or trait association. Those are rsWGCNA decisions.

## 2. Build once, reuse across projects

1. Choose a genome FASTA, GFF3, and (for reference-pack mode) protein or CDS FASTA.
2. Choose a reference pack / motif library.
3. Select mapping mode:
   - Strict for legacy-equivalent, conservative background;
   - Balanced for strict + calibrated DBD;
   - Exploratory for sequence-supported DBD candidates.
4. Set promoter window, p-value, background, threads.
5. Run once.
6. Inspect coverage, mapping statuses, testability, and target exploration.
7. Export the background manifest, TF catalogue, TF–motif table, optional site table, and TF–target prior.
8. Import the prior into one or more rsWGCNA projects.

The same background can be reused as long as the genome assembly, gene ID namespace, reference pack, promoter definition, and scan parameters are unchanged. Reuse across projects does not re-run the scan.

## 3. Strict workflow

**Use when:** you need the existing full-length RBH semantics, or a conservative prior for publication.

1. Run with `mapping_mode=strict`.
2. Review the catalogue: 271 candidate TFs, 242 testable on HF1.
3. Export prior with policy `strict`.
4. Treat all targets as level D; strict mapping does not validate the target.

## 4. Balanced workflow

**Use when:** strict leaves too many identified TFs without motifs, but you still want calibrated assignments only.

1. Run with `mapping_mode=balanced`.
2. Check that calibrated DBD rows are labeled `dbd_calibrated`, not `dbd_exploratory`.
3. Expect a small increment over strict (HF1: +9 TFs).
4. Export prior with policy `strict_calibrated`.

## 5. Exploratory workflow

**Use when:** you need candidate targets for sequence-supported TFs that strict/balanced cannot map.

1. Run with `mapping_mode=exploratory`, `sites_only=true` if you do not need a materialized prior.
2. Open the overview: check actual target coverage, ambiguous TFs, testable motifs, and background.
3. Use the TF catalogue filters to find:
   - `dbd_exploratory` TFs with targets;
   - `no_supported_motif` TFs;
   - `threshold_unattainable` motifs.
4. Open a TF profile: read the diagnostic reason, motif assignments, testability, shared-PWM warning, and candidate targets.
5. Open a target promoter: inspect strand, site sequence, and target → regulators.
6. Use enrichment only as exploratory evidence, with the displayed caveats.
7. Export the site table and TF–motif table; do not expect a TF–target prior from a site-only run.

Exploratory targets are never automatically promoted. A user who wants them in a prior must run a non-site-only exploratory job and choose policy `extended` or `selected` with `include_exploratory=true`.

## 6. Unmapped TF workflow

When a TF has no supported motif:

1. The catalogue shows `no_supported_motif` and a diagnostic reason.
2. The TF profile explains that no target panel exists.
3. The motif explorer may show a family-associated motif resource; this is not an individual TF assignment.
4. The user can switch to exploratory mode, but must not interpret a family-level motif as a validated TF–motif mapping.
5. If nothing supports the TF, it remains visible in the catalogue with zero motifs/targets; it is not silently dropped.

## 7. Target-first workflow

**Use when:** you have a gene of interest and want candidate upstream regulators.

1. Query the target or click it from a TF profile.
2. Inspect the promoter track: +/− strand sites, TSS, position, p-value.
3. Click a site for sequence and raw evidence.
4. Read the target → candidate regulators table.
5. Filter by mapping status, motif, or max p.
6. Note shared-PWM warnings; multiple TFs at one site are not independent.
7. Export evidence or use the prior in rsWGCNA for module-level context.

## 8. rsWGCNA import workflow

1. In ReguTrace, export a prior (`.db`) with a documented policy.
2. In the dashboard Data Import, upload it as `tf_target`.
3. The importer accepts the canonical SQLite database with `hits` and `pairs` tables.
4. Run the rsWGCNA steps; the prior is level D candidate evidence, not validated regulation.
5. Do not expect ReguTrace to perform expression filtering, module detection, or trait association.

## 9. Provenance workflow

Before reusing a background:

1. Open provenance.
2. Check genome/GFF/motif checksums and reference-pack version.
3. Check gene ID namespace and promoter definition.
4. Check scan parameters (p-value, background mode, exact-PWM reuse, site-only).
5. Read warnings and limitations.
6. Record the manifest with the project.

## 10. What not to do

- Do not compare target counts across mapping modes without noting that exploratory includes uncalibrated mappings.
- Do not treat `threshold_unattainable` as “no target”.
- Do not treat shared-PWM TFs as independent regulators.
- Do not use exploratory coverage as an accuracy estimate.
- Do not feed a site-only job's site table to rsWGCNA as if it were a TF–target prior.
- Do not re-run threshold optimization on the calibration dataset to make a mode look better.
