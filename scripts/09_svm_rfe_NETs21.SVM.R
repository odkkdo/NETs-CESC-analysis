# =====================================================================
# 09 - SVM-RFE feature selection (caret implementation)
# Rebuilt from the analysis working file: "Version 1/17.SVM/NETs21.SVM.R"
# REPRODUCIBILITY FIX: the original file had `#set.seed(12345)` commented out;
#   the seed is activated here so the 10-fold RFE resampling is reproducible.
#   (The separate mSVM-RFE helper NETs21.msvmRFE.R is not included: the
#   working copy contained a non-functional time-gated snippet and the
#   caret-based RFE below is the implementation that produced the published
#   feature ranking.)
# Requires: data/merge.normalize.txt, gene_lists/interGenes.txt (7 genes)
# Output: SVM-RFE.gene.txt (ranked features)
# =====================================================================

library(limma)
library(caret)
library(e1071)
set.seed(12345)          # FIXED: was commented out in the working file

outDir <- "results/09_svm_rfe"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

rt <- read.table("data/merge.normalize.txt", header = TRUE, sep = "\t", check.names = FALSE)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)

geneRT <- read.table("gene_lists/interGenes.txt", header = FALSE, sep = "\t", check.names = FALSE)
data <- data[as.vector(geneRT[, 1]), , drop = FALSE]
data <- as.data.frame(t(data))

group <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", rownames(data))
data <- cbind(group, data)
data$group <- factor(data$group, levels = c("Normal", "Cancer"))

x <- data[, -1]
y <- data$group

ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = 10)
rfeResult <- rfe(x, y, sizes = 1:ncol(x), rfeControl = ctrl)

top.features <- rfeResult$variables[, c("Variables", "Overall")]
write.table(top.features, file = file.path(outDir, "feature_svm.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

errors <- 1 - rfeResult$results$Accuracy
pdf(file = file.path(outDir, "errors.pdf"), width = 5, height = 5)
plot(1:nrow(rfeResult$results), errors, type = "b", pch = 19, xlab = "Number of Features",
     ylab = "Error Rate", main = "RFE Cross-Validation Error")
dev.off()

accuracy <- rfeResult$results$Accuracy
pdf(file = file.path(outDir, "accuracy.pdf"), width = 5, height = 5)
plot(1:nrow(rfeResult$results), accuracy, type = "b", pch = 19,
     xlab = "Variables", ylab = "Accuracy (Cross-Validation)",
     main = "Feature Selection Accuracy")
dev.off()

optimalGenes <- predictors(rfeResult)
write.table(optimalGenes, file = file.path(outDir, "SVM-RFE.gene.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
