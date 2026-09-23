args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) stop("Usage: write_r_versions.R OUTPUT.tsv")
out <- args[1]

pkgs <- c("edgeR", "data.table", "ggplot2", "patchwork", "scales")
versions <- vapply(
    pkgs,
    function(x) as.character(packageVersion(x)),
    character(1)
)

d <- data.frame(
    environment = c("r", rep("r", length(pkgs))),
    software = c("R", pkgs),
    version = c(R.version.string, versions),
    stringsAsFactors = FALSE
)

write.table(
    d,
    file=out,
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)
