# =============================================================================
# Step 0 — Check that packages and input files are in place
# Run in RStudio with the workshop folder as working directory, or:
#   Rscript scripts/00_check_setup.R
# =============================================================================
pk <- c("DESeq2", "WGCNA", "GENIE3", "igraph", "ggplot2", "ggrepel", "pheatmap")
ok <- TRUE
cat("== R packages ==\n")
for (p in pk) {
  have <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  [%s] %s\n", if (have) "ok     " else "MISSING", p))
  ok <- ok && have
}
cat("== Input files ==\n")
for (f in c("data/counts.tsv", "data/samples.tsv", "data/regulators.txt", "data/gene_annotation.tsv")) {
  have <- file.exists(f)
  cat(sprintf("  [%s] %s\n", if (have) "ok     " else "MISSING", f))
  ok <- ok && have
}
if (ok) cat("\nAll good — you are ready to start.\n") else
  cat("\nSomething is missing. Is the conda environment active, and is the working directory the workshop folder?\n",
      "Current working directory:", getwd(), "\n")
