# =====================================================================
# 05b - KEGG pathway enrichment of all significant DEGs
# Rebuilt from the analysis working file: "Version 1/6.富集/KEGG/NETs10.KEGG.R"
# NOTE: enrichment is performed on ALL 123 DEGs (up + down together).
# Requires: results/04_deg/diff.txt (script 04)
# =====================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

pvalueFilter  <- 0.05
adjPvalFilter <- 0.05

colorSel <- "p.adjust"
if (adjPvalFilter > 0.05) { colorSel <- "pvalue" }

outDir <- "results/05b_KEGG"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

rt <- read.table("results/04_deg/diff.txt", header = TRUE, sep = "\t", check.names = FALSE)

genes <- unique(as.vector(rt[, 1]))
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(entrezIDs)
rt <- cbind(rt, entrezIDs)
colnames(rt)[1] <- "gene"
rt <- rt[rt[, "entrezIDs"] != "NA", ]
gene <- rt$entrezIDs

kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1)
kk@result$Description <- gsub(" - Homo sapiens \\(human\\)", "", kk@result$Description)
KEGG <- as.data.frame(kk)
KEGG$geneID <- as.character(sapply(KEGG$geneID, function(x)
  paste(rt$gene[match(strsplit(x, "/")[[1]], as.character(rt$entrezIDs))], collapse = "/")))
KEGG <- KEGG[(KEGG$pvalue < pvalueFilter & KEGG$p.adjust < adjPvalFilter), ]
write.table(KEGG, file = file.path(outDir, "KEGG.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

showNum <- 30
if (nrow(KEGG) < showNum) { showNum <- nrow(KEGG) }

pdf(file = file.path(outDir, "barplot.pdf"), width = 8, height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, label_format = 100, color = colorSel)
dev.off()

pdf(file = file.path(outDir, "bubble.pdf"), width = 8, height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio", label_format = 100, color = colorSel)
dev.off()
