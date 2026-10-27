#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
})

args <- commandArgs(
  trailingOnly = TRUE
)

get_arg <- function(flag) {
  i <- match(
    flag,
    args
  )

  if (
    is.na(i) ||
    i == length(args)
  ) {
    stop(
      paste(
        "Missing",
        flag
      )
    )
  }

  args[i + 1]
}

split_csv <- function(x) {
  y <- strsplit(
    x,
    ",",
    fixed = TRUE
  )[[1]]

  y[nzchar(y)]
}


# ============================================================
# ARGUMENTS
# ============================================================

infile <- get_arg(
  "--input"
)

library_file <- get_arg(
  "--library-sizes"
)

outdir <- get_arg(
  "--outdir"
)

fdr_cut <- as.numeric(
  get_arg("--fdr")
)

lfc_cut <- as.numeric(
  get_arg("--min-abs-log2fc")
)

min_count <- as.numeric(
  get_arg("--min-count")
)

min_total_count <- as.numeric(
  get_arg("--min-total-count")
)

control_samples <- split_csv(
  get_arg("--control-samples")
)

test_samples <- split_csv(
  get_arg("--test-samples")
)


if (
  length(control_samples) < 2
) {
  stop(
    "At least two control samples are required"
  )
}

if (
  length(test_samples) < 2
) {
  stop(
    "At least two test samples are required"
  )
}


samples <- c(
  control_samples,
  test_samples
)

group <- factor(
  c(
    rep(
      "control",
      length(control_samples)
    ),
    rep(
      "test",
      length(test_samples)
    )
  ),
  levels = c(
    "control",
    "test"
  )
)


dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# COUNTS
# ============================================================

dt <- fread(
  infile
)

count_cols <- paste0(
  samples,
  "_raw_count"
)

if (
  !all(
    count_cols %in% names(dt)
  )
) {
  stop(
    paste(
      "Missing count columns:",
      paste(
        setdiff(
          count_cols,
          names(dt)
        ),
        collapse = ", "
      )
    )
  )
}


counts <- as.matrix(
  dt[, ..count_cols]
)

storage.mode(
  counts
) <- "double"

rownames(
  counts
) <- dt[["loop_id"]]


if (
  anyNA(counts)
) {
  stop(
    "NA values found in raw loop counts"
  )
}

if (
  any(counts < 0)
) {
  stop(
    "Negative raw loop counts found"
  )
}

if (
  any(
    counts != floor(counts)
  )
) {
  stop(
    "Non-integer raw loop counts found"
  )
}


# ============================================================
# EXTERNAL LIBRARY SIZES
# total full-depth contacts in each sample cooler
# ============================================================

lib_dt <- fread(
  library_file
)

required_lib_cols <- c(
  "sample",
  "total_contacts"
)

if (
  !all(
    required_lib_cols %in%
      names(lib_dt)
  )
) {
  stop(
    "Library-size table is malformed"
  )
}


lib_size <- lib_dt[
  match(
    samples,
    sample
  ),
  total_contacts
]


if (
  anyNA(lib_size) ||
  any(lib_size <= 0)
) {
  stop(
    "Invalid external library sizes"
  )
}


# ============================================================
# edgeR
# ============================================================

y <- DGEList(
  counts = counts,
  group = group,
  lib.size = lib_size
)


keep <- filterByExpr(
  y,
  group = group,
  min.count = min_count,
  min.total.count = min_total_count
)


filter_summary <- data.table(
  metric = c(
    "condition_blind_candidates",
    "features_retained_for_edgeR"
  ),
  value = c(
    nrow(dt),
    sum(keep)
  )
)


fwrite(
  filter_summary,
  file.path(
    outdir,
    "loop_filtering_summary.tsv"
  ),
  sep = "\t"
)


if (
  sum(keep) < 10
) {
  stop(
    "Too few features retained after edgeR filtering"
  )
}


# Keep externally supplied library sizes.
y <- y[
  keep,
  ,
  keep.lib.sizes = TRUE
]


# Composition scaling on top of total-contact exposure.
y <- calcNormFactors(
  y,
  method = "TMM"
)


normalization_dt <- data.table(
  sample = samples,
  group = as.character(group),
  library_size = y$samples$lib.size,
  TMM_norm_factor = y$samples$norm.factors,
  effective_library_size =
    y$samples$lib.size *
    y$samples$norm.factors
)


fwrite(
  normalization_dt,
  file.path(
    outdir,
    "loop_normalization_factors.tsv"
  ),
  sep = "\t"
)


design <- model.matrix(
  ~ group
)


y <- estimateDisp(
  y,
  design,
  robust = TRUE
)


fit <- glmQLFit(
  y,
  design,
  robust = TRUE
)


qlf <- glmQLFTest(
  fit,
  coef = "grouptest"
)


res <- as.data.table(
  topTags(
    qlf,
    n = Inf,
    sort.by = "none"
  )$table,
  keep.rownames = "loop_id"
)


# ============================================================
# METADATA
# ============================================================

meta_cols <- intersect(
  c(
    "loop_id",
    "chrom1",
    "start1",
    "end1",
    "chrom2",
    "start2",
    "end2",
    "distance_bp",
    "mustache_fdr",
    "mustache_detection_scale"
  ),
  names(dt)
)


meta <- dt[
  ,
  ..meta_cols
]


res <- merge(
  meta,
  res,
  by = "loop_id",
  all.y = TRUE,
  sort = FALSE
)


# ============================================================
# NORMALIZED CPM
# ============================================================

norm_cpm <- as.data.table(
  cpm(
    y,
    log = FALSE,
    normalized.lib.sizes = TRUE
  ),
  keep.rownames = "loop_id"
)


setnames(
  norm_cpm,
  c(
    "loop_id",
    samples
  )
)


res <- merge(
  res,
  norm_cpm,
  by = "loop_id",
  all.x = TRUE,
  sort = FALSE
)


res[
  ,
  control_mean_CPM :=
    rowMeans(
      .SD,
      na.rm = TRUE
    ),
  .SDcols = control_samples
]


res[
  ,
  test_mean_CPM :=
    rowMeans(
      .SD,
      na.rm = TRUE
    ),
  .SDcols = test_samples
]


# ============================================================
# CLONE CONCORDANCE
# Descriptive clone-specific effects;
# inferential test remains KO group vs WT group.
# ============================================================

pc <- 0.1

test_lfc_cols <- character()


for (
  s in test_samples
) {

  col <- paste0(
    s,
    "_log2FC_vs_control_mean"
  )

  res[
    ,
    (col) :=
      log2(
        (get(s) + pc) /
        (control_mean_CPM + pc)
      )
  ]

  test_lfc_cols <- c(
    test_lfc_cols,
    col
  )
}


same_direction <- rep(
  TRUE,
  nrow(res)
)


for (
  col in test_lfc_cols
) {

  same_direction <- (
    same_direction &
    is.finite(
      res[[col]]
    ) &
    sign(
      res[[col]]
    ) ==
    sign(
      res$logFC
    )
  )
}


res[
  ,
  all_test_same_direction :=
    same_direction
]


# ============================================================
# PRE-SPECIFIED RESULT CLASSES
# ============================================================

res[
  ,
  differential_loop :=
    FDR < fdr_cut &
    abs(logFC) >= lfc_cut
]


res[
  ,
  concordant_core :=
    differential_loop &
    all_test_same_direction
]


res[
  ,
  direction := fifelse(
    differential_loop,
    fifelse(
      logFC > 0,
      "increased_in_test",
      "decreased_in_test"
    ),
    "not_differential"
  )
]


res[
  ,
  core_direction := fifelse(
    concordant_core,
    fifelse(
      logFC > 0,
      "increased_in_test",
      "decreased_in_test"
    ),
    "not_core"
  )
]


res[
  ,
  abs_logFC := abs(logFC)
]


setorder(
  res,
  FDR,
  -abs_logFC
)


res[
  ,
  abs_logFC := NULL
]


# ============================================================
# OUTPUT TABLES
# ============================================================

fwrite(
  res,
  file.path(
    outdir,
    "differential_loops_all.tsv.gz"
  ),
  sep = "\t"
)


fwrite(
  res[
    differential_loop == TRUE
  ],
  file.path(
    outdir,
    "differential_loops_FDR_effect.tsv"
  ),
  sep = "\t"
)


fwrite(
  res[
    concordant_core == TRUE
  ],
  file.path(
    outdir,
    "differential_loops_core.tsv"
  ),
  sep = "\t"
)


summary_dt <- data.table(
  metric = c(
    "condition_blind_candidates",
    "tested_after_filter",
    "FDR_significant",
    "FDR_and_effect_size",
    "concordant_core",
    "core_increased",
    "core_decreased"
  ),
  value = c(
    nrow(dt),
    nrow(res),
    sum(
      res$FDR < fdr_cut,
      na.rm = TRUE
    ),
    sum(
      res$differential_loop,
      na.rm = TRUE
    ),
    sum(
      res$concordant_core,
      na.rm = TRUE
    ),
    sum(
      res$core_direction ==
        "increased_in_test",
      na.rm = TRUE
    ),
    sum(
      res$core_direction ==
        "decreased_in_test",
      na.rm = TRUE
    )
  )
)


fwrite(
  summary_dt,
  file.path(
    outdir,
    "differential_loop_summary.tsv"
  ),
  sep = "\t"
)


# ============================================================
# VOLCANO
# ============================================================

plot_dt <- copy(
  res
)


plot_dt[
  ,
  neglog10FDR :=
    -log10(
      pmax(
        FDR,
        1e-300
      )
    )
]


plot_dt[
  ,
  class := fifelse(
    concordant_core,
    core_direction,
    fifelse(
      differential_loop,
      "differential_nonconcordant",
      "other"
    )
  )
]


p <- ggplot(
  plot_dt,
  aes(
    x = logFC,
    y = neglog10FDR
  )
) +
  geom_point(
    aes(
      alpha = class
    ),
    size = 0.6
  ) +
  geom_hline(
    yintercept =
      -log10(fdr_cut),
    linetype = 2,
    linewidth = 0.35
  ) +
  geom_vline(
    xintercept = c(
      -lfc_cut,
      lfc_cut
    ),
    linetype = 2,
    linewidth = 0.35
  ) +
  scale_alpha_manual(
    values = c(
      "other" = 0.08,
      "differential_nonconcordant" = 0.35,
      "increased_in_test" = 0.9,
      "decreased_in_test" = 0.9
    )
  ) +
  labs(
    title = "Differential loops",
    x = "log2FC (KO / WT)",
    y = "-log10(FDR)",
    alpha = NULL
  ) +
  theme_classic(
    base_size = 10
  )


ggsave(
  file.path(
    outdir,
    "EXP_differential_loops_volcano.pdf"
  ),
  p,
  width = 6.2,
  height = 4.5
)


# ============================================================
# FIXED QC OUTPUT
# ============================================================

pdf(
  file.path(
    outdir,
    "EXP_differential_loops_QC.pdf"
  ),
  width = 7,
  height = 6
)


plotMDS(
  y,
  labels = samples,
  col = as.integer(group) + 1,
  main = "MDS - loop counts"
)


plotBCV(
  y,
  main = "Biological coefficient of variation"
)


plotMD(
  qlf,
  main = "MD plot"
)

abline(
  h = c(
    -lfc_cut,
    lfc_cut
  ),
  lty = 2
)


logcpm <- cpm(
  y,
  log = TRUE,
  prior.count = 2
)


cor_mat <- cor(
  logcpm,
  method = "pearson"
)


heatmap(
  cor_mat,
  Rowv = NA,
  Colv = NA,
  scale = "none",
  margins = c(8, 8),
  main = "Sample correlation - normalized loop counts"
)


dev.off()


print(
  summary_dt
)

print(
  normalization_dt
)
