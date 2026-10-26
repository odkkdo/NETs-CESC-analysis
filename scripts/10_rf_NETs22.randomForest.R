# =====================================================================
# 10 - Random forest feature importance
# Rebuilt from the analysis working file: "Version 1/18.RF/NETs22.randomForest.R"
# (seed 12345 was already present in the working file)
# Requires: data/merge.normalize.txt, gene_lists/interGenes.txt (7 genes)
# Output: RF.gene.txt -> importance score > 5 (top 3: GZMB, ARG1, EPCAM)
# =====================================================================

library(limma)
library(randomForest)
library(ggplot2)
set.seed(12345)

outDir <- "results/10_rf"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

rt <- read.table("data/merge.normalize.txt", header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)

geneRT <- read.table("gene_lists/interGenes.txt", header = FALSE, sep = "\t", check.names = FALSE)
data <- data[as.vector(geneRT[, 1]), , drop = FALSE]
data <- t(data)
group <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(data))

rf <- randomForest(as.factor(group) ~ ., data = data, ntree = 500)
pdf(file = file.path(outDir, "forest.pdf"), width = 5, height = 5)
plot(rf, main = "Random forest", lwd = 2)
dev.off()

optionTrees <- which.min(rf$err.rate[, 1])
rf2 <- randomForest(as.factor(group) ~ ., data = data, ntree = optionTrees)

importance <- importance(x = rf2)

rfGenes <- importance[order(importance[, "MeanDecreaseGini"], decreasing = TRUE), ]
rfGenes <- names(rfGenes[rfGenes > 5])
write.table(rfGenes, file = file.path(outDir, "RF.gene.txt"),
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

importance <- importance[order(importance[, "MeanDecreaseGini"], decreasing = TRUE), , drop = FALSE]
rt <- cbind(row.names(importance), as.data.frame(importance))
colnames(rt) <- c("ID", "Importance")
if (nrow(rt) > 30) { rt <- rt[1:30, ] }
rt$ID <- factor(rt$ID, levels = rev(rt$ID))
rt$Importance <- as.numeric(rt$Importance)
p <- ggplot(rt, aes(Importance, ID)) +
  geom_segment(aes(x = 0, xend = Importance, y = ID, yend = ID), color = "grey", cex = 1.5) +
  geom_point(aes(color = Importance), cex = 3.5) +
  scale_colour_gradient(low = "#5050FFFF", high = "#BB0021FF") +
  labs(x = "Importance", y = "") +
  theme(axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10)) +
  theme_bw()

pdf(file = file.path(outDir, "geneImportance.pdf"), width = 5, height = 4.5)
print(p)
dev.off()
