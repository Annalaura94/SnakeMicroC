#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) {
    stop(paste("Missing", flag))
  }
  args[i + 1]
}

split_csv <- function(x) {
  y <- strsplit(x, ",", fixed = TRUE)[[1]]
  y[nzchar(y)]
}

infile <- get_arg("--input")
outdir <- get_arg("--outdir")
fdr_cut <- as.numeric(get_arg("--fdr"))
lfc_cut <- as.numeric(get_arg("--min-abs-log2fc"))
min_count <- as.integer(get_arg("--min-count"))
min_total_count <- as.integer(get_arg("--min-total-count"))

control_samples <- split_csv(get_arg("--control-samples"))
test_samples <- split_csv(get_arg("--test-samples"))

if (length(control_samples) < 2) {
  stop("At least two control samples are required")
}

if (length(test_samples) < 2) {
  stop("At least two test samples are required")
}

samples <- c(control_samples, test_samples)

group <- factor(
  c(
    rep("control", length(control_samples)),
    rep("test", length(test_samples))
  ),
  levels = c("control", "test")
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(infile)

count_cols <- samples

if (!all(count_cols %in% names(dt))) {
  stop(
    paste(
      "Missing count columns:",
      paste(setdiff(count_cols, names(dt)), collapse = ", ")
    )
  )
}

counts <- as.matrix(dt[, ..count_cols])
storage.mode(counts) <- "integer"
rownames(counts) <- dt[["contact_id"]]

y <- DGEList(counts = counts, group = group)

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
  value = c(nrow(dt), sum(keep))
)

fwrite(
  filter_summary,
  file.path(outdir, "contact50_edgeR_filtering_summary.tsv"),
  sep = "\t"
)

if (sum(keep) < 100) {
  stop("Too few features retained after edgeR filtering")
}

y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

design <- model.matrix(~ group)

y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = "grouptest")

res <- as.data.table(
  topTags(
    qlf,
    n = Inf,
    sort.by = "none"
  )$table,
  keep.rownames = "contact_id"
)

meta_cols <- c(
  "contact_id",
  "chrom1", "start1", "end1",
  "chrom2", "start2", "end2",
  "distance_bp"
)

meta <- dt[, ..meta_cols]

res <- merge(
  meta,
  res,
  by = "contact_id",
  all.y = TRUE,
  sort = FALSE
)

norm_cpm <- as.data.table(
  cpm(
    y,
    log = FALSE,
    normalized.lib.sizes = TRUE
  ),
  keep.rownames = "contact_id"
)

setnames(
  norm_cpm,
  c("contact_id", samples)
)

res <- merge(
  res,
  norm_cpm,
  by = "contact_id",
  all.x = TRUE,
  sort = FALSE
)

res[, control_mean_CPM :=
  rowMeans(.SD, na.rm = TRUE),
  .SDcols = control_samples
]

res[, test_mean_CPM :=
  rowMeans(.SD, na.rm = TRUE),
  .SDcols = test_samples
]

pc <- 0.1

test_lfc_cols <- character()

for (s in test_samples) {
  col <- paste0(s, "_log2FC_vs_control_mean")

  res[, (col) :=
    log2((get(s) + pc) / (control_mean_CPM + pc))
  ]

  test_lfc_cols <- c(test_lfc_cols, col)
}

same_direction <- rep(TRUE, nrow(res))

for (col in test_lfc_cols) {
  same_direction <- (
    same_direction &
    is.finite(res[[col]]) &
    sign(res[[col]]) == sign(res$logFC)
  )
}

res[, all_test_same_direction := same_direction]

res[, differential_contact_core :=
  FDR < fdr_cut &
  abs(logFC) >= lfc_cut &
  all_test_same_direction
]

res[, direction := fifelse(
  differential_contact_core,
  fifelse(
    logFC > 0,
    "increased_in_test",
    "decreased_in_test"
  ),
  "not_core"
)]

res[, abs_logFC := abs(logFC)]
setorder(res, FDR, -abs_logFC)
res[, abs_logFC := NULL]

fwrite(
  res,
  file.path(outdir, "differential_contacts50_all.tsv.gz"),
  sep = "\t"
)

fwrite(
  res[get("differential_contact_core") == TRUE],
  file.path(outdir, "differential_contacts50_core.tsv"),
  sep = "\t"
)

summary_dt <- data.table(
  metric = c(
    "tested_after_filter",
    "FDR_significant",
    "core_differential_features",
    "core_increased",
    "core_decreased"
  ),
  value = c(
    nrow(res),
    sum(res$FDR < fdr_cut, na.rm = TRUE),
    sum(res[["differential_contact_core"]], na.rm = TRUE),
    sum(res$direction == "increased_in_test", na.rm = TRUE),
    sum(res$direction == "decreased_in_test", na.rm = TRUE)
  )
)

fwrite(
  summary_dt,
  file.path(outdir, "differential_contacts50_summary.tsv"),
  sep = "\t"
)

plot_dt <- copy(res)
plot_dt[, neglog10FDR := -log10(pmax(FDR, 1e-300))]
plot_dt[, class := fifelse(
  get("differential_contact_core"),
  direction,
  "other"
)]

p <- ggplot(
  plot_dt,
  aes(x = logFC, y = neglog10FDR)
) +
  geom_point(aes(alpha = class), size = 0.6) +
  geom_hline(
    yintercept = -log10(fdr_cut),
    linetype = 2,
    linewidth = 0.35
  ) +
  geom_vline(
    xintercept = c(-lfc_cut, lfc_cut),
    linetype = 2,
    linewidth = 0.35
  ) +
  scale_alpha_manual(
    values = c(
      "other" = 0.08,
      "increased_in_test" = 0.9,
      "decreased_in_test" = 0.9
    )
  ) +
  labs(
    title = "Differential contacts",
    x = "log2FC (test / control)",
    y = "-log10(FDR)",
    alpha = NULL
  ) +
  theme_classic(base_size = 10)

ggsave(
  file.path(outdir, "EXP_independent_contacts50_volcano.pdf"),
  p,
  width = 6.2,
  height = 4.5
)

print(summary_dt)
