# Phase 4 UI audit — ReguTrace toolkit page

**Date:** 2026-11-25
**Scope:** `wgcna-dashboard/frontend/src/views/ToolkitRegutraceView.vue` and its ReguTrace components.

## 1. Audit method

Before changing anything, the existing page was exercised with a real HF1 exploratory site-only job:

- entered from the toolkit route with `job`, `target`, and `tf` query parameters;
- checked initial load, run configuration, progress, query, TF view, promoter view, export buttons, diagnostics, and the existing TF evidence/structure panels;
- inspected the Vue component tree and API client types;
- checked that evidence was already disk-indexed rather than loaded into the browser.

The audit found that the page already had a solid visual language (Bootstrap cards, small muted descriptions, badges, compact tables, ECharts/structure tracks) and a disk-backed evidence query path. Phase 4 therefore extended the page rather than redesigning it.

## 2. What existed

- upload/reference-pack configuration, run controls, progress, cancellation;
- free query by TF or target, preview table, promoter track, raw evidence rows;
- TF structure/domain panel, sequence panel, expression-independent evidence;
- evidence/promoter/annotation exports;
- diagnostics panel with input checks and AI-assistant handoff;
- summary tiles for genes, motifs, evidence rows, etc.

## 3. Gaps found

1. **Coverage was not decomposable.** A user could see “N TFs” and “N evidence rows” but not whether a TF lacked a motif, lacked a testable motif, or lacked a significant hit.
2. **Ambiguous, untestable, and no-hit states were invisible** in the catalogue.
3. **TF profile had no per-TF target list** with p-value ordering and pagination.
4. **No target → regulator aggregation**; only raw evidence rows in the target direction.
5. **Motif registry was not explorable**; no real PWM logo and no testability detail.
6. **No enrichment surface** for the existing exploratory Fisher/BH implementation.
7. **No provenance view** exposing checksums, pack version, scan manifest, or limitations.
8. **Background/prior export policies were not exposed** in the UI.
9. **Mapping confidence could be confused with target confidence** because statuses were not shown next to targets.

## 4. What was added

### Overview

- four coverage layers: TF identification, motif mapping, testability, target prediction;
- mapping-status counts (strict, calibrated, ambiguous, exploratory, no supported motif);
- background mode and actual A/C/G/T frequencies;
- scan runtime, cache hits/misses, candidate sites, site-only label;
- actual target coverage and denominator.

### Catalogue

- filters for mapping status and testability;
- columns for motif count, target count, testable/unattainable motifs, mapping status, and diagnostic reason;
- compact status badges; annotation and domain architecture remain available.

### TF profile

- header with mapping status, method, target count, shared-PWM count;
- diagnostic reason for unmapped or untestable TFs;
- motif assignment table with testability and minimum attainable p;
- real PWM sequence logo computed from the matrix;
- paginated target table with motif/status/max-p filters;
- shared-PWM warning;
- lazy protein sequence.

### Target view

- promoter track (existing) with +/− strand lanes, TSS, site table, sequence;
- raw evidence rows (existing);
- **target → candidate regulators** table, paginated, with mapping statuses and shared-PWM warning.

### Motif explorer

- search by motif/family/donor, filter by testability;
- real PWM logo;
- PWM hash, minimum attainable p, requested cutoff, background;
- mapped TFs and specificity groups;
- registry provenance table (donor, source, method, species, license).

### Enrichment

- gene set, optional background, motif IDs, optional max p;
- Fisher/BH result table;
- explicit caveats: binary promoter hit, no GC/length matching, no gene-set overlap or motif-similarity correction, shared-PWM dependence, exploratory only.

### Provenance

- input checksums and sizes;
- reference pack key/version/species/checksums;
- scan manifest and pipeline report;
- warnings and limitations.

### Export controls

- prior policy selector (`strict`, `strict_calibrated`, `extended`, `selected`);
- background downloads: manifest, TF catalogue, TF–motif table, site table;
- site-only warning explaining that no TF–target prior is materialized.

## 5. Design-language decisions

- Reused Bootstrap cards, badges, tables, and existing spacing; no new design system.
- Kept descriptions short and muted; technical caveats are visible but not dominant.
- Used `text-bg-success/primary/warning/info/secondary` consistently for mapping status, not for biological “truth”.
- Sequence logos are generated from the actual normalized PWM matrix; no placeholder or fabricated logo.
- Long lists are paginated or truncated with a `title` tooltip; the UI never renders thousands of rows.
- Numbers use locale formatting; p-values use scientific notation below `1e-3`.
- All new asynchronous views use request sequence guards so a slow response cannot overwrite a newer selection.

## 6. Verification

- `vue-tsc -b && vite build` passes.
- Browser screenshots captured with a real HF1 job at 1440 px width:
  - [`overview-final.png`](figures/phase4/overview-final.png)
  - [`catalogue-final.png`](figures/phase4/catalogue-final.png)
  - [`tf-profile-final.png`](figures/phase4/tf-profile-final.png)
  - [`motif-detail-final.png`](figures/phase4/motif-detail-final.png)
  - [`promoter-final.png`](figures/phase4/promoter-final.png)
  - [`regulators-final.png`](figures/phase4/regulators-final.png)
  - [`enrichment-final2.png`](figures/phase4/enrichment-final2.png)
  - [`provenance-final.png`](figures/phase4/provenance-final.png)
- The promoter screenshot shows real +/− strand lanes; the regulator screenshot shows 625 candidate TFs for one target with shared-PWM warnings.

## 7. Remaining UX limitations

- The initial job summary is ~1.08 MB because it embeds the full TF catalogue. It loads once; a future version could split catalogue pagination.
- The first target aggregation for a TF with 20k–36k targets takes 2.8–6.3 s. Subsequent pages are cached.
- The legacy free-query path can return a 483 kB payload and take ~4 s for a large TF.
- No virtualized table component is used; pagination caps rows instead.
- No dedicated frontend unit tests; verification is type-check, build, API integration, and real-browser screenshots.
