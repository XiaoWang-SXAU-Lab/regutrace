# Phase 4 performance measurements

**Date:** 2026-11-25
**Environment:** Windows 11, local SSD, Chrome headless (CDP), Go server v3, Rust `regutrace` release binary, HF1 data set.

## 1. Data set

| Item | Value |
|---|---:|
| Genes / promoters | 37,415 |
| Promoter length | 1,000 bp upstream, 0 downstream |
| Motif library | JASPAR 2024, 927 motifs |
| Unique PWM groups | 677 |
| Unique PWM sites | 8,191,211 |
| Motif-ID site hits | 8,203,471 |
| Mapped TFs | 1,143 |
| TF–motif mapping rows | 32,792 |
| Site TSV size | 1.99 GB |

## 2. Scan performance

Site-only exploratory scan, 8 threads, uniform background, exact-PWM reuse, warm Phase 3 cache:

| Metric | Result |
|---|---:|
| Wall time | ~135 s |
| Unique PWM sites | 8,191,211 |
| Cache hits / misses | 215 / 397 |
| Peak working set (500 ms sampling) | 498 MB |
| Peak private memory | 580 MB |
| Evidence rows emitted | 0 (site-only) |
| Site TSV write | 1.99 GB |
| Promoter FASTA | 38 MB |
| Testability table | 168 kB |
| Promoter summary | 6.2 MB |

The dashboard job reported a scan runtime of 144.1 s; the full job (upload → scan → site index → success) took ~146–151 s, including the site index build.

## 3. Backend memory

Sampled every 2 s:

| Phase | Working set |
|---|---:|
| Server before job | ~360 MB |
| After job success, index built | ~625 MB |
| After repeated large exploratory queries | peak ~2,005 MB |
| Idle after heavy queries | ~2,005 MB (Go heap retained) |

The 2.0 GB peak is real for this workload and is the main memory cost of the site-only design. It comes from 8.2 M target offsets, 8.2 M motif offsets, per-key slice overhead, target bitsets, and transient aggregation maps for TFs with 20k–36k targets. A future version could persist a compact binary index or use memory-mapped arrays; that is not implemented in Phase 4.

**Operational note:** the dashboard keeps the job registry in memory. Restarting the server makes an existing job ID unavailable even though its files remain on disk; the job must be rerun. This was observed during performance testing.

## 4. API latency and payload

Measured with `curl` against the live job, after the index was built.

| Endpoint | Condition | Time | Payload |
|---|---|---:|---:|
| `GET /jobs/{id}` | summary, 1,981 catalogue rows | 5 ms | 1,082,657 B |
| `GET /targets` | 1,277-target TF, cold | 20 ms | 9,868 B |
| `GET /targets` | 21,808-target TF, cold | 2.80 s | 24,298 B |
| `GET /targets` | 31,669-target TF, cold | 3.58 s | 18,499 B |
| `GET /targets` | 36,513-target TF, cold | 6.26 s | 15,783 B |
| `GET /targets` | same TF, cached page | 1–2 ms | 15,783 B |
| `GET /targets` | cached offset 36,000 | 1.3 ms | 10,328 B |
| `GET /regulators` | 625 TFs for one target | 79 ms | 11,922 B |
| `GET /promoter` | 191 sites | 664 ms | 501,212 B |
| `GET /motifs` | 25 rows | 4 ms | 26,426 B |
| `GET /motifs/{id}` | one motif | 102 ms | 9,943 B |
| `GET /provenance` | warm | 2 ms | 5,197 B |
| `POST /enrich` | 2 motifs, 3 genes | 242 ms | 1,456 B |
| `GET /query` | legacy free query, large TF | 4.15 s | 483,093 B |

The bounded TF target cache changes repeated pages from seconds to milliseconds. It is keyed by TF, status filter, max p, and mapping mode; it holds at most four results.

## 5. Browser initial load

Chrome headless, 1440 px viewport, warm Vite dev server, cold median TF (21,808 targets):

| Metric | Result |
|---|---:|
| DOMContentLoaded | 71 ms |
| Load event | 80 ms |
| TF view API | 9 ms |
| Motif list API | 7 ms |
| Target API (cold) | 3,255 ms |
| Promoter API | 921 ms |
| Regulators API | 88 ms |
| JS heap used | 47.7 MB |
| JS heap total | 50.6 MB |
| DOM nodes | 4,479 |

The page shell is interactive in <100 ms. The target table is the last major block to populate, at ~3.3 s for a median-sized TF.

A second measurement with the 36,513-target TF (cached from prior requests) showed 43.1 MB heap, 3,698 nodes, and 1–3 ms target/regulator calls.

## 6. Background export sizes

| Artifact | Size | Time |
|---|---:|---:|
| TF catalogue | 317 kB | 5 ms |
| TF–motif (exploratory) | 6.9 MB | 34 ms |
| TF–motif (balanced) | 47.7 kB | 3 ms |
| TF–motif (strict) | 44.8 kB | 4 ms |
| Manifest (catalogue excluded) | 3,660 B | 479 ms |
| Site table | 1.99 GB | 3.4 s stream |

The manifest originally embedded the full 1.7 MB catalogue; it now carries slim statistics and points to the separate catalogue file.

## 7. Interpretation

- Scan time and memory are acceptable for a 1.0–1.5 GB plant genome on a workstation.
- The site index is the memory bottleneck: ~2 GB peak under heavy exploration is high but avoids a much larger row materialization.
- Cold target aggregation is the main interaction bottleneck. It is linear in the number of sites for the TF's motifs, not in genome size alone.
- Caching makes pagination and repeated inspection fast without unbounded memory growth.
- The legacy free-query path remains slower than the dedicated TF/target endpoints; the UI uses the dedicated endpoints for profiles.

## 8. Not measured

- Cold OS page cache / first read after reboot.
- Network-attached storage or spinning disk.
- Multiple concurrent users.
- Production (minified, non-Vite) bundle timing.
- Full rsWGCNA run with a Phase 4 prior.
