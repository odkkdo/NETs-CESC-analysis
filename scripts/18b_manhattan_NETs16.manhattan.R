# =====================================================================
# 18b - Manhattan plot of DEG association P values
# Rebuilt from the analysis working file: "Version 1/12.曼哈顿图/NETs16.manhattan.R"
# Requires: results/04_deg/all.txt, reference_data/geneREF.txt,
#           gene_lists/interGenes.txt
# =====================================================================

library(dplyr)
library(tidyr)
library(CMplot)

outDir <- "results/18b_manhattan"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

posRT <- read.table("reference_data/geneREF.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
diffRT <- read.table("results/04_deg/all.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

sameGene <- intersect(row.names(posRT), row.names(diffRT))
rt <- cbind(posRT[sameGene, ], diffRT[sameGene, ])

data <- rt[, c("Chr", "Start", "P.Value")]
data <- cbind(Gene = row.names(data), data)
colnames(data) <- c("Gene", "CHR", "BP", "pvalue")
data$CHR <- gsub("chr", "", data$CHR)

geneRT <- read.table("gene_lists/interGenes.txt", header = FALSE, sep = "\t", check.names = FALSE)
interGene <- intersect(sameGene, as.vector(geneRT[, 1]))

CMplot(data, plot.type = "m", LOG10 = TRUE, chr.den.col = NULL,
       highlight = interGene, highlight.text = interGene,
       highlight.cex = 1.5, highlight.text.cex = 1.2,
       file = "pdf", file.output = TRUE, width = 8, height = 5, verbose = TRUE)
