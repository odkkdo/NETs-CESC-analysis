# =====================================================================
# 15 - Consensus clustering of the 40 tumor samples
# Rebuilt from the analysis working file: "Version 1/24.cluster/Cluster.R"
# Parameters used in the study: Euclidean distance + kmeans algorithm,
#   reps = 50, pItem = 0.8, seed = 123456, K evaluated 2-12 (as reported in
#   the manuscript; the surviving working file had maxK = 9 - set to 12 here
#   to reproduce the reported range). K = 3 was selected ->
#   C1 = 20, C2 = 12, C3 = 8 (matches the shipped cluster assignments).
# Requires: data/merge.normalize.txt, gene_lists/interFeatureGenes.txt
# =====================================================================

library(limma)
library(ConsensusClusterPlus)

outDir <- "results/15_cluster"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

rt <- read.table("data/merge.normalize.txt", header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)

geneRT <- read.table("gene_lists/interFeatureGenes.txt", header = FALSE, sep = "\t", check.names = FALSE)
sameGene <- intersect(as.vector(geneRT[, 1]), row.names(data))
data <- data[sameGene, ]

# ---- keep tumor samples only ----
Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data <- data[, Type == "Cancer"]

# ---- consensus clustering ----
maxK <- 12
results <- ConsensusClusterPlus(data,
                                maxK = maxK,
                                reps = 50,
                                pItem = 0.8,
                                pFeature = 1,
                                title = outDir,
                                clusterAlg = "km",
                                distance = "euclidean",
                                seed = 123456,
                                plot = "png")

# ---- export the K = 3 solution ----
clusterNum <- 3                       # K = 3 selected in the study
cluster <- results[[clusterNum]][["consensusClass"]]
cluster <- as.data.frame(cluster)
colnames(cluster) <- c("Cluster")
cluster$Cluster <- paste0("C", cluster$Cluster)
outTab <- cbind(t(data), cluster)
outTab <- rbind(ID = colnames(outTab), outTab)
write.table(outTab, file = file.path(outDir, "cluster.txt"), sep = "\t", quote = FALSE, col.names = FALSE)
