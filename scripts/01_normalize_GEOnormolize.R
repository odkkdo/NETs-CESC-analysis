# =====================================================================
# 01 - Per-dataset normalization (limma normalizeBetweenArrays)
# Rebuilt from the analysis working file: "Version 1/2.处理样本名/GEOnormolize.R"
#
# Input  (per dataset, in a dedicated folder):
#   geneMatrix.txt : rows = genes (col 1 = gene symbol), cols 2+ = samples
#   s1.txt         : sample names of the control/normal group (one per line)
#   s2.txt         : sample names of the tumor group (one per line)
# Output:
#   <geoID>.normalize.txt  (column names "<GSM-id>_<Type>")
#
# NOTE: run once per dataset. In this study:
#   GSE63514  (platform GPL570)   -> kept 24 Normal + 28 Cancer
#   GSE192897 (platform GPL27956, NanoString PanCancer IO 360)
#                                 -> kept 14 Normal + 12 Cancer
#   (sample-retention rules are described in the manuscript and README)
# =====================================================================

library(limma)

geoID     <- "GSE63514"        # change to "GSE192897" for the second dataset
inputFile <- "geneMatrix.txt"
conFile   <- "s1.txt"
treatFile <- "s2.txt"

# ---- read expression matrix ----
rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
rt <- avereps(data)            # average duplicated gene symbols

# ---- log2-transform if the matrix is not already log-scale ----
qx   <- as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- ((qx[5] > 100) || ((qx[6] - qx[1]) > 50 && qx[2] > 0))
if (LogC) {
  rt[rt < 0] <- 0
  rt <- log2(rt + 1)
}
data <- normalizeBetweenArrays(rt)

# ---- subset to the retained samples ----
sample1 <- read.table(conFile,   header = FALSE, sep = "\t", check.names = FALSE)
sample2 <- read.table(treatFile, header = FALSE, sep = "\t", check.names = FALSE)
sampleName1 <- gsub("^ | $", "", as.vector(sample1[, 1]))
sampleName2 <- gsub("^ | $", "", as.vector(sample2[, 1]))
conData   <- data[, sampleName1]
treatData <- data[, sampleName2]
data      <- cbind(conData, treatData)
conNum    <- ncol(conData)
treatNum  <- ncol(treatData)

Type    <- c(rep("Normal", conNum), rep("Cancer", treatNum))
outData <- rbind(id = paste0(colnames(data), "_", Type), data)
write.table(outData, file = paste0(geoID, ".normalize.txt"),
            sep = "\t", quote = FALSE, col.names = FALSE)
