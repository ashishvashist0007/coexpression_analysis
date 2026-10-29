# =============================================================================
# Step 5 — Guilt by association: predict what unknown genes do
#
# "Tell me who your friends are, and I'll tell you who you are."
# A gene whose expression moves with many cell-cycle genes is probably
# involved in the cell cycle too. About 1,500 of our genes have no known role.
#
# Method A — module vote:  an unknown gene inherits the most common known role
#                          in its WGCNA module.
# Method B — friend vote:  an unknown gene's 10 most correlated genes ("friends")
#                          vote with their known roles.
# Then we test method B honestly: hide some genes whose role we DO know, predict
# them, and measure how often we get them right (cross-validation, AUROC).
#
# In:   results/01_prepare, results/02_modules, data/gene_annotation.tsv
# Out:  results/05_gba/  (predictions, accuracy, a "friends" plot)
# =============================================================================
suppressPackageStartupMessages({ library(WGCNA); library(igraph); library(ggplot2) })
options(stringsAsFactors = FALSE)
out <- "results/05_gba"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

prep <- readRDS("results/01_prepare/prepared.rds")
nw   <- readRDS("results/02_modules/network.rds")
datExpr <- prep$datExpr
genes  <- colnames(datExpr)
module <- setNames(nw$moduleColors, genes)
role   <- read.delim("data/gene_annotation.tsv", row.names = 1)[genes, "known_role"]
names(role) <- genes
table(role)


## Method A — module vote -----------------------------------------------------------------
mod_role <- sapply(split(role, module), function(r) {
  r <- r[r != "unknown"]
  if (length(r) < 3) return(NA)
  tab <- sort(table(r), decreasing = TRUE)
  if (tab[1] / length(r) < 0.6) NA else names(tab)[1]    # need a clear majority
})
mod_role["grey"] <- NA                                   # grey = no module, no vote
mod_role


## Method B — vote of the 10 closest "friends" -------------------------------------------------
K <- 10
C <- WGCNA::cor(datExpr)                                 # gene x gene correlation
diag(C) <- NA
friends <- t(apply(C, 1, function(x) names(sort(x, decreasing = TRUE))[1:K]))

friend_vote <- function(known) {                         # known: named role vector ("unknown" = no vote)
  roles <- setdiff(unique(known), "unknown")
  votes <- t(apply(friends, 1, function(f) table(factor(known[f], levels = roles))))
  votes / K                                              # fraction of friends with each role
}
V <- friend_vote(role)

best <- colnames(V)[max.col(V, ties.method = "first")]
pred <- data.frame(gene = genes, module = module,
                   module_vote = unname(mod_role[module]),
                   friend_vote = ifelse(apply(V, 1, max) >= 0.5, best, NA),
                   friend_score = apply(V, 1, max),
                   top_friends = apply(friends[, 1:5], 1, paste, collapse = ","))
unknown_pred <- subset(pred[role == "unknown", ], !is.na(friend_vote) | !is.na(module_vote))
unknown_pred <- unknown_pred[order(-unknown_pred$friend_score), ]
cat("Unknown genes with a prediction:", nrow(unknown_pred), "\n")
head(unknown_pred, 20)
table(friend_vote = unknown_pred$friend_vote, useNA = "ifany")
write.table(unknown_pred, file.path(out, "predictions_unknown_genes.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)


## Does it work? Cross-validation on genes we DO know -------------------------------------------
# For each role: hide one third of its known genes, predict them from the rest,
# and check whether the hidden ones get higher scores than everyone else.
# AUROC: 0.5 = no better than a coin flip, 1.0 = perfect.
auroc <- function(score, positive) {
  r <- rank(score)
  n1 <- sum(positive); n0 <- sum(!positive)
  (sum(r[positive]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
set.seed(42)
cv <- do.call(rbind, lapply(setdiff(unique(role), "unknown"), function(r) {
  pos <- names(role)[role == r]
  fold <- sample(rep(1:3, length.out = length(pos)))
  a <- sapply(1:3, function(f) {
    hidden <- pos[fold == f]
    known <- role; known[hidden] <- "unknown"
    s <- friend_vote(known)[, r]
    test <- names(role)[!(names(role) %in% setdiff(pos, hidden))]   # hidden positives + all non-members
    auroc(s[test], test %in% hidden)
  })
  data.frame(role = r, known_genes = length(pos), AUROC = round(mean(a), 3))
}))
cv[order(-cv$AUROC), ]
write.table(cv, file.path(out, "cross_validation_auroc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


## Picture: one unknown gene and its friends --------------------------------------------------------
example <- "MKI67"      # the classic proliferation marker -- "unknown" in our annotation table
if (role[example] != "unknown")                          # (in case you changed the annotation)
  example <- unknown_pred$gene[1]
f <- friends[example, ]
g <- make_star(K + 1, mode = "undirected")
V(g)$name  <- c(example, f)
pal <- c("cell cycle" = "#2F6FA8", "interferon response" = "#C2621B", "hypoxia response" = "#7A5BA6",
         "cholesterol synthesis" = "#3E8E5A", "heat shock response" = "#B5443B", "ribosome" = "#C9A227",
         "unknown" = "#B8C0C6")
V(g)$color <- pal[c("unknown", role[f])]
E(g)$width <- 1 + 8 * (C[example, f] - min(C[example, f])) / diff(range(C[example, f]))
png(file.path(out, paste0("friends_of_", example, ".png")), width = 1000, height = 900, res = 130)
set.seed(1)
plot(g, layout = layout_as_star(g), vertex.size = c(30, rep(20, K)), vertex.label.cex = 0.75,
     vertex.label.color = "#13222B", vertex.frame.color = NA, edge.color = "grey60",
     main = paste0(example, ": role unknown. Its 10 closest friends vote ",
                   sum(role[f] == "cell cycle"), "/10 for 'cell cycle'"))
legend("bottomleft", legend = names(pal), fill = pal, bty = "n", cex = 0.75)
invisible(dev.off())

saveRDS(list(pred = pred, V = V, cv = cv, friends = friends), file.path(out, "gba.rds"))
