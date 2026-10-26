# =====================================================================
# 02 - Merge the two normalized datasets and remove batch effects (ComBat)
# Rebuilt from the analysis working file: "Version 1/3.整合 merge.normalize/geo06.sva.R"
#
# BUG FIX vs the original working file: inside the merge loop the skip
#   conditions referenced `file` (the loop variable of the *gene-list* loop
#   above) instead of `inputFile`, so the guard `if(file=="merge.preNorm.txt"){next}`
#   never matched during the merge pass. Corrected below by using an explicit
#   input list instead of a directory scan.
#
# Input : GSE63514.normalize.txt / GSE192897.normalize.txt (from script 01)
# Output: merge.preNorm.txt    (merged, before batch correction)
#         merge.normalize.txt  (merged, ComBat-corrected; 78 samples = 38 Normal + 40 Cancer)
# The final merge.normalize.txt is shipped with this repository in data/.
# =====================================================================

library(limma)
library(sva)

inputs <- c("GSE63514.normalize.txt", "GSE192897.normalize.txt")

# ---- gene intersection across datasets ----
geneList <- list()
for (file in inputs) {
  rt <- read.table(file, header = TRUE, sep = "\t", check.names = FALSE)
  geneNames <- as.vector(rt[, 1])
  geneList[[gsub("\\.normalize\\.txt$", "", file)]] <- unique(geneNames)
}
interGenes <- Reduce(intersect, geneList)

# ---- merge per-dataset matrices on the shared genes ----
allTab <- data.frame(); batchType <- c()
for (i in seq_along(inputs)) {
  inputFile <- inputs[i]
  rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE)
  rt <- as.matrix(rt)
  rownames(rt) <- rt[, 1]
  exp <- rt[, 2:ncol(rt)]
  dimnames <- list(rownames(exp), colnames(exp))
  data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
  rt <- avereps(data)
  colnames(rt) <- paste0(gsub("\\.normalize\\.txt$", "", inputFile), "_", colnames(rt))

  if (i == 1) {
    allTab <- rt[interGenes, ]
  } else {
    allTab <- cbind(allTab, rt[interGenes, ])
  }
  batchType <- c(batchType, rep(i, ncol(rt)))
}

outTab <- rbind(geneNames = colnames(allTab), allTab)
write.table(outTab, file = "merge.preNorm.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# ---- ComBat batch correction ----
outTab <- ComBat(allTab, batchType, par.prior = TRUE)
outTab <- rbind(geneNames = colnames(outTab), outTab)
write.table(outTab, file = "merge.normalize.txt", sep = "\t", quote = FALSE, col.names = FALSE)
