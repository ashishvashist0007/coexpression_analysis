# =============================================================================
# Step 1 — From counts to a clean expression matrix
#
#  1. load the raw counts and the sample traits
#  2. drop genes with too few reads
#  3. variance-stabilise (DESeq2) so every gene is on a comparable log-like scale
#  4. cluster the samples and remove outliers
#
# In:   data/counts.tsv, data/samples.tsv
# Out:  results/01_prepare/  (sample tree plot, prepared data for steps 2-5)
#
# Run section by section in RStudio (working directory = workshop folder),
# or all at once:  Rscript scripts/01_prepare_data.R
# =============================================================================
suppressPackageStartupMessages({ library(DESeq2); library(WGCNA) })
options(stringsAsFactors = FALSE)
out <- "results/01_prepare"
dir.create(out, recursive = TRUE, showWarnings = FALSE)


## 1. Load ----------------------------------------------------------------------
counts <- as.matrix(read.delim("data/counts.tsv", row.names = 1, check.names = FALSE))
traits <- read.delim("data/samples.tsv", row.names = 1)
stopifnot(identical(colnames(counts), rownames(traits)))
dim(counts)          # genes x samples
head(traits)


## 2. Filter genes with too few reads ------------------------------------------------
# Lowly expressed genes are mostly noise, and noise makes random correlations.
keep <- rowSums(counts >= 10) >= 0.9 * ncol(counts)
table(keep)
counts <- counts[keep, ]


## 3. Variance-stabilising transformation -------------------------------------------
# Correlation needs values on a comparable scale: raw counts are dominated by
# highly expressed genes, and their variance grows with the mean.
# (For a full transcriptome use the faster vst(); with ~1,600 genes we use the
#  exact function it approximates.)
dds <- DESeqDataSetFromMatrix(counts, colData = traits, design = ~ 1)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "local")
expr <- assay(vsd)                    # genes x samples, log2-like scale

# Real data: you would usually keep the 5,000-10,000 most variable genes here.
# Our dataset is small, so we keep every expressed gene.

# WGCNA wants samples in rows and genes in columns
datExpr <- t(expr)
gsg <- goodSamplesGenes(datExpr, verbose = 0)
gsg$allOK


## 4. Look for outlier samples ---------------------------------------------------------
sampleTree <- hclust(dist(datExpr), method = "average")
trait_cols <- numbers2colors(data.frame(oxygen = traits$oxygen_pct, ifn = traits$ifn_treated,
                                        growth = traits$growth_rate,
                                        batch = as.numeric(traits$batch == "B")), signed = FALSE)

cut_height <- 40          # chosen by looking at the tree: change it if your tree differs
png(file.path(out, "sample_tree.png"), width = 1100, height = 700, res = 120)
plotDendroAndColors(sampleTree, trait_cols, groupLabels = c("oxygen", "IFN", "growth", "batch"),
                    main = "Sample clustering: is any sample an outlier?", cex.dendroLabels = 0.8)
abline(h = cut_height, col = "red")
invisible(dev.off())

clust <- cutreeStatic(sampleTree, cutHeight = cut_height, minSize = 5)
table(clust)                                      # cluster 0 = samples above the red line
outliers <- rownames(datExpr)[clust != 1]
cat("Outlier samples removed:", if (length(outliers)) outliers else "none", "\n")

datExpr <- datExpr[clust == 1, ]
traits  <- traits[rownames(datExpr), ]
cat("Kept", nrow(datExpr), "samples and", ncol(datExpr), "genes\n")

saveRDS(list(datExpr = datExpr, traits = traits, outliers = outliers),
        file.path(out, "prepared.rds"))
