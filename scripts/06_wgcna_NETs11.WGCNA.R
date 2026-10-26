# =====================================================================
# 06 - WGCNA co-expression network
# Rebuilt from the analysis working file: "Version 1/7.WGCNA/NETs11.WGCNA.R"
# Parameters used in the study: SD > 0.5 filter (634 genes retained); no
#   outlier samples removed (cutHeight 20000 kept all 78); unsigned network;
#   soft threshold beta = 8 (first power reaching scale-free R^2 > 0.9,
#   R^2 ~ 0.92); minModuleSize = 50; module merge height 0.30; traits =
#   Normal / Cancer 0-1 indicator columns.
# MEturquoise (148 genes) correlated most strongly with the Cancer trait
#   (r = 0.37, P = 8e-4) and was carried forward; MEblue had |r| = 0.78 but
#   pointed towards the Normal trait (both reported transparently in the
#   response letter).
# =====================================================================

library(limma)
library(WGCNA)

expFile <- "data/merge.normalize.txt"
outDir <- "results/06_WGCNA"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

# ---- read and filter ----
rt <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)
data <- data[apply(data, 1, sd) > 0.5, ]     # keep genes with SD > 0.5

Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data <- data[, order(Type)]
Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
conCount   <- length(Type[Type == "Normal"])
treatCount <- length(Type[Type == "Cancer"])
datExpr0 <- t(data)

# ---- check missing values ----
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# ---- sample clustering / outlier check ----
sampleTree <- hclust(dist(datExpr0), method = "average")
pdf(file = file.path(outDir, "1_sample_cluster.pdf"), width = 18, height = 10)
par(cex = 0.6); par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 20000, col = "red")
dev.off()

clust <- cutreeStatic(sampleTree, cutHeight = 20000, minSize = 10)
table(clust)
keepSamples <- (clust == 1)
datExpr0 <- datExpr0[keepSamples, ]

# ---- trait data (0/1 indicators) ----
traitData <- data.frame(Normal = c(rep(1, conCount), rep(0, treatCount)),
                        Cancer = c(rep(0, conCount), rep(1, treatCount)))
row.names(traitData) <- colnames(data)
fpkmSamples <- rownames(datExpr0)
traitSamples <- rownames(traitData)
sameSample <- intersect(fpkmSamples, traitSamples)
datExpr0 <- datExpr0[sameSample, ]
datTraits <- traitData[sameSample, ]

sampleTree2 <- hclust(dist(datExpr0), method = "average")
traitColors <- numbers2colors(datTraits, signed = FALSE)
pdf(file = file.path(outDir, "2_sample_heatmap.pdf"), width = 18, height = 12)
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# ---- soft threshold ----
enableWGCNAThreads()
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
pdf(file = file.path(outDir, "3_scale_independence.pdf"), width = 9, height = 5)
par(mfrow = c(1, 2)); cex1 <- 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

softPower <- sft$powerEstimate      # beta = 8 in the original run (first R^2 > 0.9)
adjacency <- adjacency(datExpr0, power = softPower)

# ---- TOM and gene clustering ----
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")
pdf(file = file.path(outDir, "4_gene_clustering.pdf"), width = 12, height = 9)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# ---- dynamic module detection ----
minModuleSize <- 50
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2,
                             pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)
pdf(file = file.path(outDir, "5_Dynamic_Tree.pdf"), width = 8, height = 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# ---- merge similar modules (height 0.30) ----
MEList <- moduleEigengenes(datExpr0, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
pdf(file = file.path(outDir, "6_Clustering_module.pdf"), width = 7, height = 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
MEDissThres <- 0.30
abline(h = MEDissThres, col = "red")
dev.off()

merge <- mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs
pdf(file = file.path(outDir, "7_merged_dynamic.pdf"), width = 9, height = 6)
plotDendroAndColors(geneTree, mergedColors, "Merged dynamic", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
moduleColors <- mergedColors
table(moduleColors)
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder) - 1
MEs <- mergedMEs

# ---- module-trait relationships ----
nGenes <- ncol(datExpr0)
nSamples <- nrow(datExpr0)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
pdf(file = file.path(outDir, "8_Module_trait.pdf"), width = 6, height = 5.5)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits), yLabels = names(MEs),
               ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.75,
               zlim = c(-1, 1), main = "Module-trait relationships")
dev.off()

y <- datTraits[, 1]
GS1 <- as.numeric(cor(y, datExpr0, use = "p"))
GeneSignificance <- abs(GS1)
ModuleSignificance <- tapply(GeneSignificance, mergedColors, mean, na.rm = TRUE)
pdf(file = file.path(outDir, "9_GeneSignificance.pdf"), width = 11, height = 7)
plotModuleSignificance(GeneSignificance, mergedColors)
dev.off()

# ---- export module gene lists ----
for (mod in 1:nrow(table(moduleColors))) {
  modules <- names(table(moduleColors))[mod]
  probes <- colnames(datExpr0)
  inModule <- (moduleColors == modules)
  modGenes <- probes[inModule]
  write.table(modGenes, file = file.path(outDir, paste0("module_", modules, ".txt")),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
