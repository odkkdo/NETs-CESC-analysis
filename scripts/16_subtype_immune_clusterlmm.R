# =====================================================================
# 16 - Immune-cell signature differences between molecular subtypes
# Rebuilt from the analysis working file: "Version 1/27.clusterlmm/clusterlmm.R"
# (no methodological changes)
# Requires: results/15_cluster/cluster.txt, results/07_ssGSEA/ssGSEA.result.txt
# =====================================================================

library(limma)
library(reshape2)
library(ggpubr)

outDir <- "results/16_subtype_immune"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

immune <- read.table("results/07_ssGSEA/ssGSEA.result.txt", header = TRUE, sep = "\t",
                     check.names = FALSE, row.names = 1)

# ---- tumor samples only ----
immune <- t(immune)
group <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(immune))
data <- immune[group == "Cancer", , drop = FALSE]

Cluster <- read.table("results/15_cluster/cluster.txt", header = TRUE, sep = "\t",
                      check.names = FALSE, row.names = 1)

sameSample <- intersect(row.names(data), row.names(Cluster))
rt <- cbind(data[sameSample, , drop = FALSE], Cluster[sameSample, "Cluster", drop = FALSE])
rt <- rt[order(rt$Cluster, decreasing = FALSE), ]

# ---- boxplot across subtypes ----
plotData <- melt(as.data.frame(rt), id.vars = c("Cluster"))
colnames(plotData) <- c("Cluster", "Immune", "Expression")

group <- levels(factor(plotData$Cluster))
bioCol <- c("#0066FF", "#FF0000", "#FF9900", "#6E568C", "#7CC767", "#223D6C",
            "#D20A13", "#FFD121", "#088247", "#11AA4D")
bioCol <- bioCol[1:length(group)]

boxplot <- ggboxplot(plotData, x = "Immune", y = "Expression", color = "Cluster",
                     xlab = "", ylab = "Fraction", legend.title = "Cluster",
                     add = "point", width = 0.8, palette = bioCol) +
  rotate_x_text(50) +
  stat_compare_means(aes(group = Cluster),
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "")),
                     label = "p.signif")

pdf(file = file.path(outDir, "immune.diff.pdf"), width = 8, height = 6)
print(boxplot)
dev.off()
