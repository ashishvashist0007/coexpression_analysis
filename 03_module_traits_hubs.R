# =============================================================================
# Step 3 — What do the modules mean?
#
#  1. module eigengenes vs sample traits (the module-trait heatmap)
#  2. what is in each module (known gene roles)
#  3. hub genes: the most connected genes of each module (kME)
#  4. gene significance vs module membership for one module
#  5. export one module for Cytoscape
#
# In:   results/01_prepare/prepared.rds, results/02_modules/network.rds,
#       data/gene_annotation.tsv
# Out:  results/03_traits_hubs/
# =============================================================================
suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
out <- "results/03_traits_hubs"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

prep <- readRDS("results/01_prepare/prepared.rds")
nw   <- readRDS("results/02_modules/network.rds")
datExpr <- prep$datExpr; traits <- prep$traits; moduleColors <- nw$moduleColors
annot <- read.delim("data/gene_annotation.tsv", row.names = 1)


## 1. Module eigengenes vs traits ------------------------------------------------------
# The eigengene is the first principal component of a module: one "average"
# expression profile that summarises all its genes.
MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
MEs <- MEs[, colnames(MEs) != "MEgrey", drop = FALSE]
traitNum <- data.frame(oxygen_pct = traits$oxygen_pct, ifn_treated = traits$ifn_treated,
                       growth_rate = traits$growth_rate, batch_B = as.numeric(traits$batch == "B"),
                       row.names = rownames(traits))
modTraitCor <- cor(MEs, traitNum, use = "p")
modTraitP   <- corPvalueStudent(modTraitCor, nrow(datExpr))
round(modTraitCor, 2)

txt <- paste0(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")")
dim(txt) <- dim(modTraitCor)
png(file.path(out, "module_trait_heatmap.png"), width = 1000, height = 900, res = 120)
par(mar = c(6, 9, 3, 2))
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitNum), yLabels = colnames(MEs),
               ySymbols = colnames(MEs), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = txt, setStdMargins = FALSE, cex.text = 0.75, zlim = c(-1, 1),
               main = "Module-trait relationships (r and p-value)")
invisible(dev.off())


## 2. What is inside each module? ---------------------------------------------------------
# Count the KNOWN roles among each module's genes ("unknown" not shown).
role <- annot[colnames(datExpr), "known_role"]
comp <- table(module = moduleColors, role = role)
comp <- comp[, colnames(comp) != "unknown", drop = FALSE]
comp
sizes <- table(moduleColors)
summary_tab <- data.frame(module = names(sizes), genes = as.integer(sizes),
  top_known_role = apply(comp[names(sizes), , drop = FALSE], 1,
                         function(x) if (sum(x) == 0) "-" else names(which.max(x))),
  annotated_with_it = apply(comp[names(sizes), , drop = FALSE], 1, max))
summary_tab


## 3. Hub genes -----------------------------------------------------------------------------
# kME = correlation of a gene with its module eigengene ("module membership").
# The genes with the highest kME are the module's hubs.
kME <- signedKME(datExpr, MEs)
colnames(kME) <- sub("^kME", "", colnames(kME))
hubs <- do.call(rbind, lapply(setdiff(unique(moduleColors), "grey"), function(m) {
  g <- colnames(datExpr)[moduleColors == m]
  k <- kME[g, m]
  o <- order(-k)[1:min(5, length(g))]
  data.frame(module = m, gene = g[o], kME = round(k[o], 3), known_role = annot[g[o], "known_role"])
}))
hubs
write.table(hubs, file.path(out, "hub_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


## 4. Gene significance vs module membership -------------------------------------------------
# For the module most correlated with growth: do the genes most central to the
# module (high kME) also track growth best (high gene significance)?
best <- sub("^ME", "", rownames(modTraitCor)[which.max(abs(modTraitCor[, "growth_rate"]))])
GS <- abs(cor(datExpr, traitNum$growth_rate, use = "p"))[, 1]
inMod <- moduleColors == best
png(file.path(out, "GS_vs_MM.png"), width = 800, height = 700, res = 120)
verboseScatterplot(abs(kME[inMod, best]), GS[inMod], col = best,
                   xlab = paste("Module membership (kME) in", best, "module"),
                   ylab = "Gene significance for growth rate",
                   main = paste("Module membership vs gene significance\n"), abline = TRUE)
invisible(dev.off())


## 5. Export one module for Cytoscape ---------------------------------------------------------
g <- colnames(datExpr)[inMod]
TOM <- TOMsimilarityFromExpr(datExpr[, g], power = nw$beta, networkType = "signed", TOMType = "signed",
                             verbose = 0)
dimnames(TOM) <- list(g, g)
cyt <- exportNetworkToCytoscape(TOM, threshold = quantile(TOM[upper.tri(TOM)], 0.9),
                                nodeNames = g, nodeAttr = role[inMod],
                                edgeFile = file.path(out, paste0("cytoscape_edges_", best, ".txt")),
                                nodeFile = file.path(out, paste0("cytoscape_nodes_", best, ".txt")))
cat("Exported the", best, "module:", nrow(cyt$edgeData), "edges between", length(g), "genes\n")

saveRDS(list(MEs = MEs, kME = kME, modTraitCor = modTraitCor, summary = summary_tab),
        file.path(out, "traits_hubs.rds"))
write.table(summary_tab, file.path(out, "module_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
