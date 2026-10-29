# =============================================================================
# Step 2 — Build the WGCNA network and find modules
#
#  1. choose the soft-thresholding power (beta)
#  2. build the network and cut it into modules (one call: blockwiseModules)
#  3. look at the dendrogram and the module sizes
#
# In:   results/01_prepare/prepared.rds
# Out:  results/02_modules/  (soft-threshold plot, dendrogram, module table)
# =============================================================================
suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
allowWGCNAThreads(2)
out <- "results/02_modules"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

prep <- readRDS("results/01_prepare/prepared.rds")
datExpr <- prep$datExpr


## 1. Pick the soft-thresholding power --------------------------------------------
# adjacency = ((1 + correlation) / 2) ^ beta     ("signed" network)
# A higher beta keeps strong correlations and shrinks weak ones towards zero.
# We look for the lowest beta where the network is roughly scale-free
# (fit index R^2 above ~0.8) while genes still keep some connections.
powers <- c(1:10, seq(12, 20, 2))
invisible(capture.output(
  sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)))
fit <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
print(data.frame(power = powers, R2 = round(fit, 3),
                 mean_connectivity = round(sft$fitIndices$mean.k., 1)))

beta <- powers[which(fit >= 0.8)[1]]
if (is.na(beta)) beta <- 12                  # fallback if the fit never reaches 0.8
cat("Chosen soft-thresholding power:", beta, "\n")

png(file.path(out, "soft_threshold.png"), width = 1200, height = 520, res = 120)
par(mfrow = c(1, 2))
plot(powers, fit, type = "n", xlab = "Soft threshold (power)", ylab = "Scale-free fit (R²)",
     main = "Scale independence", ylim = c(min(fit, 0), 1))
text(powers, fit, labels = powers, col = "#C2621B")
abline(h = 0.8, col = "#2F6FA8", lty = 2)
plot(powers, sft$fitIndices$mean.k., type = "n", xlab = "Soft threshold (power)",
     ylab = "Mean connectivity", main = "Mean connectivity")
text(powers, sft$fitIndices$mean.k., labels = powers, col = "#C2621B")
invisible(dev.off())


## 2. Build the network and detect modules -----------------------------------------
# blockwiseModules does it all: adjacency -> topological overlap (TOM) ->
# hierarchical clustering -> dynamic tree cut -> merge similar modules.
net <- blockwiseModules(datExpr, power = beta,
                        networkType = "signed", TOMType = "signed",
                        minModuleSize = 10,         # smallest module allowed
                        mergeCutHeight = 0.25,      # merge modules whose eigengenes correlate > 0.75
                        maxBlockSize = 5000, numericLabels = FALSE,
                        saveTOMs = FALSE,
                        randomSeed = 1234, verbose = 0)
moduleColors <- net$colors
table(moduleColors)                           # "grey" = genes in no module


## 3. Plot the gene dendrogram with module colours ------------------------------------
png(file.path(out, "gene_dendrogram_modules.png"), width = 1300, height = 650, res = 120)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], "Module",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and modules")
invisible(dev.off())

modules <- data.frame(gene = colnames(datExpr), module = moduleColors)
write.table(modules, file.path(out, "gene_modules.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(net = net, beta = beta, moduleColors = moduleColors), file.path(out, "network.rds"))
