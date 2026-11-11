#!/usr/bin/env Rscript
# Fair benchmark: regutrace vs Biostrings::matchPWM.
#
# A single relative threshold (e.g. 80% of max score) is not comparable across
# implementations. This script calibrates thresholds on negative promoters and
# reports precision/recall at a matched false-positive rate (FPR), plus average precision (AP).
#
# Synthetic design: n/2 positive promoters each with one planted consensus site
# (plus or minus strand), n/2 negative promoters with no planted site.
#
# Usage:
#   Rscript scripts/benchmark.R <regutrace> <jaspar.txt> <motif_id> <outdir> [n] [len] [real_promoters.fa]
#
# Example:
#   Rscript scripts/benchmark.R target/release/regutrace.exe \
#     _rsWGCNA-local/regutrace-data/JASPAR2024_CORE_plants_non-redundant_pfms_jaspar.txt \
#     MA1758.2 _rsWGCNA-local/regutrace-benchmark 1200 1000 \
#     _rsWGCNA-local/regutrace-benchmark/real_promoters.fa

suppressPackageStartupMessages(library(Biostrings))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: benchmark.R <regutrace> <jaspar> <motif_id> <outdir> [n] [len] [real_promoters.fa]")
}
bin <- normalizePath(args[1], mustWork = TRUE)
jaspar <- normalizePath(args[2], mustWork = TRUE)
motif_id <- args[3]
outdir <- args[4]
n <- if (length(args) >= 5) as.integer(args[5]) else 1200L
promoter_len <- if (length(args) >= 6) as.integer(args[6]) else 1000L
real_fa <- if (length(args) >= 7) args[7] else NA_character_
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

bg <- c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)
target_fpr <- 0.01

read_pfm <- function(path, id) {
  lines <- readLines(path, warn = FALSE)
  i <- grep(paste0("^>", id, "([[:space:]]|$)"), lines)
  if (length(i) != 1) stop("motif not found or duplicated: ", id)
  rows <- list()
  k <- i + 1L
  while (length(rows) < 4L && k <= length(lines)) {
    line <- trimws(lines[k])
    if (nzchar(line) && substr(line, 1, 1) %in% c("A", "C", "G", "T")) {
      b <- substr(line, 1, 1)
      clean <- gsub("[", "", substr(line, 2, nchar(line)), fixed = TRUE)
      clean <- gsub("]", "", clean, fixed = TRUE)
      rows[[b]] <- scan(text = clean, quiet = TRUE)
    }
    k <- k + 1L
  }
  m <- do.call(rbind, rows[c("A", "C", "G", "T")])
  storage.mode(m) <- "integer"
  m
}

pfm <- read_pfm(jaspar, motif_id)
p_our <- (pfm + 0.5) / (colSums(pfm) + 2)
lo_our <- log2(p_our / bg)
max_our <- sum(apply(lo_our, 2, max))
consensus <- paste(apply(p_our, 2, function(z) c("A", "C", "G", "T")[which.max(z)]),
                   collapse = "")
motif_len <- nchar(consensus)
rc <- function(s) paste(rev(chartr("ACGT", "TGCA", strsplit(s, "")[[1]])), collapse = "")

# Biostrings needs equal PFM column sums; add the tiny deficit to A.
pfm_eq <- pfm
target <- max(colSums(pfm_eq))
for (j in seq_len(ncol(pfm_eq))) {
  pfm_eq["A", j] <- pfm_eq["A", j] + (target - sum(pfm_eq[, j]))
}
pwm <- PWM(matrix(as.integer(pfm_eq), nrow = 4,
                  dimnames = list(c("A", "C", "G", "T"), NULL)),
           type = "log2probratio", prior.params = bg)

set.seed(42)
rand_seq <- function(n) {
  paste(sample(c("A", "C", "G", "T"), n, replace = TRUE, prob = bg), collapse = "")
}
n_pos <- n %/% 2L
n_neg <- n - n_pos
n_consensus <- n_pos %/% 2L
n_weak <- n_pos - n_consensus
promoters <- character(n)
ids <- paste0(c(rep("pos_", n_pos), rep("neg_", n_neg)), seq_len(n))
truth <- data.frame(promoter = character(), start = integer(), end = integer(),
                    strand = character(), stringsAsFactors = FALSE)

score_site <- function(site) {
  bases <- strsplit(site, "")[[1]]
  idx <- match(bases, c("A", "C", "G", "T"))
  sum(lo_our[cbind(idx, seq_along(idx))])
}
# Weak sites: sampled from the PWM, constrained to 50-85% of the maximum score.
weak_site <- function() {
  for (k in seq_len(500)) {
    site <- paste(apply(p_our, 2, function(prob) sample(c("A", "C", "G", "T"), 1, prob = prob)),
                  collapse = "")
    s <- score_site(site)
    if (s >= 0.5 * max_our && s <= 0.85 * max_our) return(site)
  }
  consensus
}

for (i in seq_len(n)) {
  seq_chars <- strsplit(rand_seq(promoter_len), "")[[1]]
  if (i <= n_pos) {
    strand <- if (runif(1) < 0.5) "+" else "-"
    site <- if (i <= n_consensus) consensus else weak_site()
    if (strand == "-") site <- rc(site)
    start <- sample.int(promoter_len - motif_len + 1L, 1)
    seq_chars[start:(start + motif_len - 1L)] <- strsplit(site, "")[[1]]
    truth <- rbind(truth, data.frame(promoter = ids[i], start = start,
                                     end = start + motif_len - 1L, strand = strand))
  }
  promoters[i] <- paste(seq_chars, collapse = "")
}
names(promoters) <- ids
fasta <- file.path(outdir, "synthetic_promoters.fa")
con <- file(fasta, "w")
for (i in seq_along(promoters)) {
  writeLines(paste0(">", ids[i]), con)
  writeLines(promoters[i], con)
}
close(con)
write.table(truth, file.path(outdir, "synthetic_truth.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

tf_map <- file.path(outdir, "benchmark_tf_map.tsv")
writeLines(c("tf_id\ttf_symbol\tfamily\tmotif_id",
             paste("TF1", "TF1", "TEST", motif_id, sep = "\t")), tf_map)

# Run regutrace at a permissive score cutoff; sweep thresholds in R.
low_our <- max_our * 0.5
our_file <- file.path(outdir, "our_hits.tsv")
run_our <- function(out_file, min_score) {
  status <- system2(bin, c("scan", "--promoters", fasta, "--motifs", jaspar,
                           "--tf-map", tf_map, "--pvalue", "1",
                           "--min-score", format(min_score, digits = 12),
                           "--out", out_file), stdout = FALSE, stderr = FALSE)
  if (!identical(status, 0L)) stop("regutrace scan failed with status ", status)
}
t_our <- replicate(3, system.time(run_our(our_file, low_our))[["elapsed"]])
our <- read.delim(our_file, check.names = FALSE)
our_hits <- data.frame(promoter = our$target_id, start = our$hit_start,
                       end = our$hit_end, strand = our$strand,
                       score = our$score, stringsAsFactors = FALSE)

# Biostrings scan with scores; also sweep thresholds in R.
scan_biostrings <- function(promoters, min_score = "50%") {
  hits <- list()
  for (i in seq_along(promoters)) {
    s <- DNAString(promoters[i])
    plus <- matchPWM(pwm, s, min.score = min_score, with.score = TRUE)
    if (length(plus) > 0) {
      hits[[length(hits) + 1L]] <- data.frame(
        promoter = names(promoters)[i], start = start(plus), end = end(plus),
        strand = "+", score = mcols(plus)$score, stringsAsFactors = FALSE)
    }
    rcs <- reverseComplement(s)
    minus <- matchPWM(pwm, rcs, min.score = min_score, with.score = TRUE)
    if (length(minus) > 0) {
      n <- nchar(promoters[i])
      hits[[length(hits) + 1L]] <- data.frame(
        promoter = names(promoters)[i],
        start = n - start(minus) - width(minus) + 2L,
        end = n - start(minus) + 1L,
        strand = "-", score = mcols(minus)$score, stringsAsFactors = FALSE)
    }
  }
  if (length(hits) == 0) return(data.frame(promoter = character(), start = integer(),
                                            end = integer(), strand = character(),
                                            score = numeric()))
  do.call(rbind, hits)
}
scan_bi <- function() scan_biostrings(promoters, "50%")
t_bi <- replicate(3, system.time(scan_bi())[["elapsed"]])
bi_hits <- scan_bi()

key <- function(d) paste(d$promoter, d$start, d$end, d$strand, sep = "|")
eval_at <- function(hits, truth, n_neg, threshold) {
  h <- hits[hits$score >= threshold, , drop = FALSE]
  hk <- key(h); tk <- key(truth)
  tp <- length(intersect(hk, tk))
  fp <- length(setdiff(hk, tk))
  fn <- length(setdiff(tk, hk))
  data.frame(threshold = threshold, TP = tp, FP = fp, FN = fn,
             FPR = fp / n_neg,
             precision = if (tp + fp > 0) tp / (tp + fp) else NA_real_,
             recall = if (tp + fn > 0) tp / (tp + fn) else NA_real_,
             F1 = if (tp > 0 && (tp + fp) > 0 && (tp + fn) > 0)
               2 * tp / (2 * tp + fp + fn) else NA_real_)
}

sweep <- function(hits, thresholds, n_neg) {
  do.call(rbind, lapply(thresholds, function(t) eval_at(hits, truth, n_neg, t)))
}
our_curve <- sweep(our_hits, seq(0.5, 1.0, length.out = 101) * max_our, n_neg)
bi_curve <- sweep(bi_hits, seq(0.5, 1.0, length.out = 101), n_neg)

choose_at_fpr <- function(curve) {
  ok <- curve[!is.na(curve$FPR) & curve$FPR <= target_fpr, , drop = FALSE]
  if (nrow(ok) == 0) ok <- curve[which.min(curve$FPR), , drop = FALSE]
  ok[which.max(ok$recall), , drop = FALSE]
}
our_sel <- choose_at_fpr(our_curve)
bi_sel <- choose_at_fpr(bi_curve)

average_precision <- function(curve) {
  x <- curve[!is.na(curve$recall) & !is.na(curve$precision), , drop = FALSE]
  if (nrow(x) < 2) return(NA_real_)
  x <- x[order(x$threshold, decreasing = TRUE), , drop = FALSE]
  ap <- 0
  prev <- 0
  for (i in seq_len(nrow(x))) {
    r <- x$recall[i]
    if (r > prev) {
      ap <- ap + (r - prev) * x$precision[i]
      prev <- r
    }
  }
  ap
}

summary <- data.frame(
  tool = c("regutrace", "Biostrings::matchPWM"),
  selected_threshold = c(our_sel$threshold, bi_sel$threshold),
  threshold_scale = c("raw log2 odds", "0-1 unit scale"),
  TP = c(our_sel$TP, bi_sel$TP),
  FP = c(our_sel$FP, bi_sel$FP),
  FN = c(our_sel$FN, bi_sel$FN),
  FPR = c(our_sel$FPR, bi_sel$FPR),
  precision = c(our_sel$precision, bi_sel$precision),
  recall = c(our_sel$recall, bi_sel$recall),
  F1 = c(our_sel$F1, bi_sel$F1),
  average_precision = c(average_precision(our_curve), average_precision(bi_curve)),
  n_positive = n_pos,
  n_negative = n_neg,
  n_consensus = n_consensus,
  n_weak = n_weak,
  target_FPR = target_fpr,
  median_seconds = c(median(t_our), median(t_bi)),
  stringsAsFactors = FALSE
)
write.table(summary, file.path(outdir, "benchmark_synthetic.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(rbind(cbind(tool = "regutrace", our_curve),
                  cbind(tool = "Biostrings::matchPWM", bi_curve)),
            file.path(outdir, "benchmark_curves.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# Real promoter concordance at each tool's matched-FPR threshold.
real_row <- NULL
if (!is.na(real_fa) && file.exists(real_fa)) {
  real <- readDNAStringSet(real_fa)
  real_names <- names(real)
  real_proms <- setNames(as.character(real), real_names)
  real_file <- file.path(outdir, "real_hits.tsv")
  run_real <- function() {
    status <- system2(bin, c("scan", "--promoters", real_fa, "--motifs", jaspar,
                             "--tf-map", tf_map, "--pvalue", "1",
                             "--min-score", format(our_sel$threshold, digits = 12),
                             "--out", real_file), stdout = FALSE, stderr = FALSE)
    if (!identical(status, 0L)) stop("regutrace scan failed with status ", status)
  }
  t_real_our <- replicate(3, system.time(run_real())[["elapsed"]])
  rh <- read.delim(real_file, check.names = FALSE)
  rh <- data.frame(promoter = rh$target_id, start = rh$hit_start, end = rh$hit_end,
                   strand = rh$strand, stringsAsFactors = FALSE)
  scan_real_bi <- function() scan_biostrings(real_proms, bi_sel$threshold)
  t_real_bi <- replicate(3, system.time(scan_real_bi())[["elapsed"]])
  rh_bi <- scan_real_bi()
  j <- length(intersect(key(rh), key(rh_bi))) / max(1, length(union(key(rh), key(rh_bi))))
  real_row <- data.frame(real_promoters = length(real), our_hits = nrow(rh),
                         biostrings_hits = nrow(rh_bi), jaccard = j,
                         our_median_seconds = median(t_real_our),
                         biostrings_median_seconds = median(t_real_bi))
  write.table(real_row, file.path(outdir, "benchmark_real.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
}

cat("\n=== matched-FPR synthetic benchmark ===\n")
print(summary, row.names = FALSE)
if (!is.null(real_row)) {
  cat("\n=== real promoter concordance ===\n")
  print(real_row, row.names = FALSE)
}
cat("\nResults:", outdir, "\n")
