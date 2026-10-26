# =====================================================================
# 08 - LASSO logistic regression feature selection
# Rebuilt from the analysis working file: "Version 1/16.lasso/NETs20.lasso.R"
# REPRODUCIBILITY FIX: the original file had `#set.seed(12345)` commented out;
#   the seed is activated here so the 10-fold cross-validation in cv.glmnet is
#   reproducible. Cosmetic repair: pd0000f( -> pdf(.
# Requires: data/merge.normalize.txt, gene_lists/interGenes.txt (7 genes)
# Output: LASSO.gene.txt -> the 3 signature genes (ARG1, GZMB, EPCAM)
# =====================================================================

library(limma)
library(glmnet)
set.seed(12345)          # FIXED: was commented out in the working file

outDir <- "results/08_lasso"
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
rt <- t(data)

x <- as.matrix(rt)
y <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(rt))
fit <- glmnet(x, y, family = "binomial", alpha = 1)

pdf(file = file.path(outDir, "lasso.pdf"), width = 6, height = 5.5)
plot(fit)
dev.off()

cvfit <- cv.glmnet(x, y, family = "binomial", alpha = 1, type.measure = "deviance", nfolds = 10)
pdf(file = file.path(outDir, "cvfit.pdf"), width = 6, height = 5.5)   # FIXED: was "pd0000f("
plot(cvfit)
dev.off()

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
lassoGene <- row.names(coef)[index]
lassoGene <- lassoGene[-1]     # drop the intercept
write.table(lassoGene, file = file.path(outDir, "LASSO.gene.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
