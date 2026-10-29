# =============================================================================
# Step 4 — Who regulates whom? Regulatory network inference with GENIE3
#
# WGCNA tells us which genes move TOGETHER (undirected).
# GENIE3 asks, for every target gene: which transcription factors best PREDICT
# its expression? It fits a random forest per target, using only the
# candidate regulators as inputs, and scores each regulator's importance.
# High importance = a candidate regulator -> target link (directed).
#
# In:   results/01_prepare/prepared.rds, results/02_modules/network.rds,
#       data/regulators.txt
# Out:  results/04_genie3/  (ranked links, TF summary, network plot)
# =============================================================================
suppressPackageStartupMessages({ library(GENIE3); library(igraph) })
options(stringsAsFactors = FALSE)
out <- "results/04_genie3"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

prep <- readRDS("results/01_prepare/prepared.rds")
nw   <- readRDS("results/02_modules/network.rds")
expr <- t(prep$datExpr)                               # GENIE3 wants genes x samples
module <- setNames(nw$moduleColors, colnames(prep$datExpr))


## 1. Candidate regulators and targets ---------------------------------------------------
regulators <- intersect(readLines("data/regulators.txt"), rownames(expr))
length(regulators)                                    # TFs we let GENIE3 use
# To keep it fast we use the genes that sit in a module, plus the TFs.
# (On a real dataset you would use a few thousand variable genes and a TF list
#  such as the human TF catalogue; it then runs for hours, not seconds.)
targets <- union(names(module)[module != "grey"], regulators)
length(targets)


## 2. Run GENIE3 ----------------------------------------------------------------------------
set.seed(123)                                         # random forests are random
weights <- GENIE3(expr[targets, ], regulators = regulators, nTrees = 500, nCores = 2)
links <- getLinkList(weights)                         # all regulator -> target pairs, ranked
links$regulatoryGene <- as.character(links$regulatoryGene)
links$targetGene <- as.character(links$targetGene)
links$target_module <- module[links$targetGene]
head(links, 15)
write.table(links, file.path(out, "genie3_links_ranked.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


## 3. Which TFs dominate the top of the list? -------------------------------------------------
top <- head(links, 300)
tf_summary <- as.data.frame.matrix(table(regulator = top$regulatoryGene, module = top$target_module))
tf_summary$total <- rowSums(tf_summary)
tf_summary <- tf_summary[order(-tf_summary$total), ]
tf_summary                                            # how many top-300 links per TF, by target module
write.table(cbind(regulator = rownames(tf_summary), tf_summary), file.path(out, "tf_summary_top300.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# For each module: the best regulator and how much better it is than the runner-up
best_reg <- do.call(rbind, lapply(setdiff(unique(module[targets]), "grey"), function(m) {
  g <- setdiff(intersect(targets, names(module)[module == m]), regulators)
  w <- sort(colMeans(t(weights[, g, drop = FALSE])), decreasing = TRUE)
  data.frame(module = m, genes = length(g), best_regulator = names(w)[1], mean_weight = round(w[1], 3),
             runner_up = names(w)[2], runner_up_weight = round(w[2], 3))
}))
best_reg[order(-best_reg$mean_weight), ]
write.table(best_reg, file.path(out, "best_regulator_per_module.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


## 4. Draw the top of the network -------------------------------------------------------------
edges <- head(links, 150)
g <- graph_from_data_frame(edges[, c("regulatoryGene", "targetGene")], directed = TRUE)
is_tf <- V(g)$name %in% regulators
mod_col <- module[V(g)$name]; mod_col[is.na(mod_col)] <- "grey"
V(g)$color <- ifelse(is_tf, "#13222B", mod_col)
V(g)$shape <- ifelse(is_tf, "square", "circle")
V(g)$size  <- ifelse(is_tf, 9, 4)
V(g)$label <- ifelse(is_tf, V(g)$name, NA)
set.seed(3)
png(file.path(out, "genie3_network_top150.png"), width = 1300, height = 1100, res = 130)
plot(g, layout = layout_with_fr(g), edge.arrow.size = 0.25, edge.color = "grey70",
     vertex.label.color = "#13222B", vertex.label.cex = 0.8, vertex.label.dist = 1.2,
     vertex.frame.color = NA, main = "GENIE3: top 150 regulator -> target links\n(squares = TFs, circles = targets coloured by WGCNA module)")
invisible(dev.off())

saveRDS(list(weights = weights, links = links), file.path(out, "genie3.rds"))
