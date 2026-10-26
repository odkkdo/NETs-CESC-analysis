# =====================================================================
# 19a - Expression boxplots (7 intersection genes) + gene-gene correlation
# Rebuilt from the analysis working file: "Version 1/10.boxplot/NETs14.boxplot.R"
# REPAIR of corrupted lines (the working file had been damaged by manual
#   editing; restored by obvious intent):
#   `p=g` / `pgboxplot(`            -> `p1 = ggboxplot(...)` + `print(p1)`
#   `dv.off()`                      -> `dev.off()`
#   `data=geneExp[,Type=="TreCancerop=F]` -> `data <- geneExp[, Type == "Cancer", drop=FALSE]`
#   duplicated `corrplot(cor1` + stray `.off()` -> single well-formed call
# The correlation block uses tumor samples only (Pearson, 40 tumors).
# Requires: data/merge.normalize.txt, gene_lists/interGenes.txt
# =====================================================================

library(limma)
library(corrplot)
library(circlize)
library(reshape2)
library(ggpubr)

outDir <- "results/19a_boxplot"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

rt <- read.table("data/merge.normalize.txt", header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)

geneRT <- read.table("gene_lists/interGenes.txt", header = FALSE, sep = "\t", check.names = FALSE)
geneExp <- data[as.vector(geneRT[, 1]), , drop = FALSE]

# ---- boxplot: Normal vs Cancer for each of the 7 genes ----
Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
expdf <- as.data.frame(t(geneExp))
expdf <- cbind(expdf, Type = Type)
plotData <- melt(expdf, id.vars = c("Type"))
colnames(plotData) <- c("Type", "Gene", "Expression")

p1 <- ggboxplot(plotData, x = "Gene", y = "Expression", color = "Type", fill = "Type",
                xlab = "", ylab = "Gene expression", legend.title = "",
                palette = c("#E64B35", "#4DBBD5"), add = "jitter",
                add.params = list(size = 1.6, alpha = 0.75), width = 0.65,
                outlier.shape = NA, alpha = 0.15) +
  stat_compare_means(aes(group = Type), method = "wilcox.test", label = "p.signif", size = 5) +
  rotate_x_text(45) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(face = "bold", size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.position = "top", legend.text = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.margin = margin(5.5, 10, 5.5, 5.5))
pdf(file = file.path(outDir, "boxplot.pdf"), width = 7, height = 5)
print(p1)
dev.off()

# ---- gene-gene correlation heatmap (tumor samples only) ----
corrData <- geneExp[, Type == "Cancer", drop = FALSE]   # REPAIRED line
rtc <- t(corrData)
cor1 <- cor(rtc)

pdf(file = file.path(outDir, "corrplot.pdf"), width = 7, height = 7)
corrplot(cor1,
         method = "square",
         order = "hclust",
         type = "upper",
         col = colorRampPalette(c("#4DBBD5", "white", "#E64B46"))(200),
         tl.col = "black",
         tl.cex = 0.9,
         tl.srt = 45,
         cl.cex = 0.9,
         cl.ratio = 0.2,
         addgrid.col = "grey70",
         diag = FALSE,
         mar = c(0, 0, 1, 0),
         addCoef.col = "black",
         number.cex = 0.7)
dev.off()
