# =====================================================================
# 20 - Drug enrichment (DSigDB via clusterProfiler::enricher)
# Rebuilt from the analysis working file: "Version 1/30.DrugEnrich/NETs34.drugEnrich.R"
# Input: gene_lists/interFeatureGenes.txt (3 signature genes) and a DSigDB
#   drug-target table. The DSigDB file (~50 MB, from the Enrichr/DSigDB
#   resource) is NOT redistributed here - download it and save the first two
#   columns (Term, Gene) as data/DSigDB_All_detailed.txt; see README.
# Result in the study: 14 candidate drugs (p < 0.05 and BH adj.P < 0.05),
#   including carboplatin and topotecan (EPCAM targets).
# =====================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(tidytable)
library(ggplot2)

pvalueFilter  <- 0.05
adjPvalFilter <- 0.05

drugFile <- "data/DSigDB_All_detailed.txt"   # user-provided (see header)
hubFile  <- "gene_lists/interFeatureGenes.txt"
outDir <- "results/20_drug"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

rt <- read.table(hubFile, header = FALSE, sep = "\t", check.names = FALSE)
genes <- unique(as.vector(rt[, 1]))

drugRT <- read.table(drugFile, header = TRUE, sep = "\t", check.names = FALSE,
                     quote = "", comment.char = "")
drugRT <- drugRT[, 1:2]

kk <- enricher(genes, pvalueCutoff = 1, qvalueCutoff = 1,
               minGSSize = 10, maxGSSize = 500,
               TERM2GENE = drugRT)

DRUG <- as.data.frame(kk)
DRUG <- DRUG[(DRUG$pvalue < pvalueFilter & DRUG$p.adjust < adjPvalFilter), ]
write.table(DRUG[, -c(3, 4)], file = file.path(outDir, "DRUG.enrich.xls"),
            sep = "\t", quote = FALSE, row.names = FALSE)

showNum <- ifelse(nrow(DRUG) > 30, 30, nrow(DRUG))
showDRUG <- DRUG[1:showNum, ]
showDRUG$Drug <- factor(showDRUG$Description, levels = rev(showDRUG$Description))

bar <- ggplot(data = showDRUG, aes(x = Count, y = Drug, fill = p.adjust)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_fill_distiller(palette = "Spectral", direction = 1) +
  labs(x = "Gene count", y = "", title = "Enriched drug") +
  theme_bw() + theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
                     axis.title = element_text(size = 11), axis.text = element_text(size = 10),
                     legend.title = element_text(size = 11), legend.text = element_text(size = 10))
pdf(file = file.path(outDir, "barplot.pdf"), width = 6, height = 5.5)
print(bar)
dev.off()

showDRUG$GeneRatio <- sapply(showDRUG$GeneRatio, function(x) eval(parse(text = x)))
showDRUG <- showDRUG[order(showDRUG$GeneRatio), ]
showDRUG$Drug <- factor(showDRUG$Drug, levels = showDRUG$Drug)
bubble <- ggplot(data = showDRUG, aes(x = GeneRatio, y = Drug)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_distiller(palette = "Spectral", direction = 1) +
  labs(x = "Gene ratio", y = "", title = "Enriched drug", size = "Count") +
  theme_bw() + theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
                     axis.title = element_text(size = 11), axis.text = element_text(size = 10),
                     legend.title = element_text(size = 11), legend.text = element_text(size = 10))
pdf(file = file.path(outDir, "bubble.pdf"), width = 6, height = 5.5)
print(bubble)
dev.off()

# ---- drug-gene network tables ----
DRUGn <- DRUG[, c("ID", "geneID")]
DRUG2 <- separate_rows(as_tidytable(DRUGn), geneID, sep = "/")
DRUG2 <- DRUG2 %>% group_by(geneID) %>% slice_head(n = 10)
networkTab <- cbind(as.data.frame(DRUG2), "Drug")
colnames(networkTab) <- c("Node1", "Node2", "Type")
write.table(networkTab, file = file.path(outDir, "net.network.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

drugNode <- cbind(unique(networkTab[, 1]), "Drug")
colnames(drugNode) <- c("Node", "Type")
geneNode <- cbind(unique(networkTab[, 2]), "Gene")
colnames(geneNode) <- c("Node", "Type")
nodeTab <- rbind(drugNode, geneNode)
write.table(nodeTab, file = file.path(outDir, "net.node.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
