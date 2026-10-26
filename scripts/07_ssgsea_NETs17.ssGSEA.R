# =====================================================================
# 07 - ssGSEA immune-cell signature scoring (28-signature panel)
# Rebuilt from the analysis working file: "Version 1/13.ssGSEA/NETs17.ssGSEA.R"
# In the merged gene panel, 22 of the 28 signatures met the minSize = 3
#   requirement and were scoreable; scores are min-max normalized.
# Requires: data/merge.normalize.txt, reference_data/immune.gmt
# =====================================================================

library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

expFile  <- "data/merge.normalize.txt"
gmtFile  <- "reference_data/immune.gmt"
outDir <- "results/07_ssGSEA"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

rt <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)

geneSets <- getGmt(gmtFile, geneIdType = SymbolIdentifier())

ssgseaPar <- ssgseaParam(data, geneSets, minSize = 3, maxSize = 500)
ssgseaScore <- gsva(ssgseaPar)

normalize <- function(x) { (x - min(x)) / (max(x) - min(x)) }
ssgseaScore <- normalize(ssgseaScore)

ssgseaOut <- rbind(id = colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut, file = file.path(outDir, "ssGSEA.result.txt"),
            sep = "\t", quote = FALSE, col.names = FALSE)
