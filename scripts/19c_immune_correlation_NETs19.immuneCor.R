# =====================================================================
# 19c - Gene x immune-signature correlation (Spearman, tumor samples only)
# Rebuilt from the analysis working file: "Version 1/15.immuneCor/NETs19.immuneCor.R"
# REPAIR of a corrupted line: `treat=grepl("_TreaCancerameSample, ignore.case=T)`
#   -> `treat <- grepl("_Cancer", sameSample, ignore.case = TRUE)`.
# SCOPE NOTE (as stated in the response letter): correlations are computed
#   within the 40 tumor samples only - this matches the original analysis
#   (the grepl filter selects "_Cancer" samples). 22 genes x 28 signatures
#   here; the 22 scoreable signatures x 7 genes = 154 tests are tabulated in
#   Supplementary Table S3B. A table export was added for transparency.
# Requires: data/merge.normalize.txt, gene_lists/interGenes.txt,
#           results/07_ssGSEA/ssGSEA.result.txt
# =====================================================================

library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)

outDir <- "results/19c_immuneCor"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

immune <- read.table("results/07_ssGSEA/ssGSEA.result.txt", header = TRUE, sep = "\t",
                     check.names = FALSE, row.names = 1)
immune <- t(immune)

rt <- read.table("data/merge.normalize.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

geneRT <- read.table("gene_lists/interGenes.txt", header = FALSE, sep = "\t", check.names = FALSE)
exp <- t(rt[as.vector(geneRT[, 1]), ])

# ---- restrict to tumor samples (REPAIRED filter) ----
sameSample <- intersect(row.names(immune), row.names(exp))
treat <- grepl("_Cancer", sameSample, ignore.case = TRUE)   # REPAIRED line
sameSample <- sameSample[treat]
immune <- immune[sameSample, , drop = FALSE]
exp <- exp[sameSample, , drop = FALSE]

# ---- Spearman correlations ----
outTab <- data.frame()
for (cell in colnames(immune)) {
  for (gene in colnames(exp)) {
    x <- as.numeric(immune[, cell])
    y <- as.numeric(exp[, gene])
    corT <- cor.test(x, y, method = "spearman")
    cor <- corT$estimate
    pvalue <- corT$p.value
    text <- ifelse(pvalue < 0.001, "***", ifelse(pvalue < 0.01, "**", ifelse(pvalue < 0.05, "*", "")))
    outTab <- rbind(outTab, cbind(Gene = gene, Immune = cell, cor, text, pvalue))
  }
}

outTab$cor <- as.numeric(outTab$cor)
outTab$pvalue <- as.numeric(outTab$pvalue)
write.table(outTab, file = file.path(outDir, "immuneCor.full.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)   # table export added

pdf(file = file.path(outDir, "cor.pdf"), width = 6, height = 6.5)
ggplot(outTab, aes(Gene, Immune)) +
  geom_tile(aes(fill = cor), colour = "grey", size = 1) +
  scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") +
  geom_text(aes(label = text), col = "black", size = 3) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold")) +
  labs(fill = paste0("***  p<0.001", "\n", "**  p<0.01", "\n", " *  p<0.05", "\n", "\n", "Correlation")) +
  scale_x_discrete(position = "bottom")
dev.off()
