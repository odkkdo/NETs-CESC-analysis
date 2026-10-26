# =====================================================================
# 14 - External validation (TCGA-CESC + GTEx, 404 samples)
# Rebuilt from the analysis working file: "Version 1/23.testDiff/testDiff.R"
#
# METHOD FIX (documented in the response letter, applied to the revision):
#   The original script refit the logistic model *inside* the validation set
#   (glm(y ~ ., data=rt)) and called predict(logit, newx=rt) (a typo - `newx`
#   is not a predict.glm argument; the call effectively returned the fitted
#   values). The resulting AUC = 0.994 (95% CI 0.986-1.000) is therefore an
#   in-cohort refitted, optimistically biased estimate. The corrected script
#   provides:
#   (1) PRIMARY - locked-coefficient extrapolation: the logistic model
#       coefficients trained on the 78-sample cohort are applied unchanged to
#       the external expression matrix
#       (intercept 0.8263, ARG1 -1.0028, GZMB 0.3379, EPCAM 0.1328)
#       -> AUC = 0.968 (95% CI 0.940-0.989);
#       Youden cutoff 0.812 -> Sensitivity 96.7%, Specificity 96.9%.
#   (2) SENSITIVITY - the original within-cohort refit, retained but clearly
#       labelled: AUC = 0.994 (95% CI 0.986-1.000);
#       Youden cutoff 0.629 -> Sensitivity 99.3%, Specificity 96.9%.
#   Single-gene external AUCs (pROC auto direction): ARG1 0.819 (raw
#   orientation 0.181, ARG1 is lower in tumors - pROC flips automatically),
#   GZMB 0.959, EPCAM 0.981.
#
# Input: data/GTExTCGA.merge.txt - NOT shipped (162 MB, exceeds GitHub limits).
#   Build it from UCSC Xena (TCGA-CESC tumors + adjacent normals, GTEx
#   normal cervix), rows = genes, columns named "<source>_<sample>_<Type>"
#   with Type in {Cancer, Normal}; see README "External validation data".
# =====================================================================

library(limma)
library(reshape2)
library(ggpubr)
library(glmnet)
library(pROC)

expFile <- "data/GTExTCGA.merge.txt"
geneFile <- "gene_lists/interFeatureGenes.txt"
outDir <- "results/14_external_validation"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

# ---- read external matrix ----
rt <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)

geneRT <- read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)
sameGene <- intersect(as.vector(geneRT[, 1]), row.names(data))
data <- t(data[sameGene, ])

Type <- gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt <- cbind(as.data.frame(data), Type)

# ---- expression boxplot ----
plotData <- melt(rt, id.vars = c("Type"))
colnames(plotData) <- c("Type", "Gene", "Expression")

p1 <- ggboxplot(plotData, x = "Gene", y = "Expression", color = "Type", fill = "Type",
                xlab = "", ylab = "Gene expression", legend.title = "",
                palette = c("#4DBBD5", "#E64B35"), add = "jitter",
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
pdf(file = file.path(outDir, "boxplot.pdf"), width = 6, height = 4.5)
print(p1)
dev.off()

# ---- per-gene ROC ----
bioCol <- c("#F05C3BFF", "#5C88DAFF", "#5CB85CFF", "#EEA236FF", "#9632B8FF", "#17BECFFF", "#BCBD22FF")
bioCol <- bioCol[1:nrow(geneRT)]
aucText <- c()
k <- 0
y <- ifelse(rt$Type == "Normal", 0, 1)
for (x in colnames(rt)[1:(ncol(rt) - 1)]) {
  k <- k + 1
  roc1 <- roc(y, as.numeric(rt[, x]))    # pROC auto-reverses the direction
  if (k == 1) {
    pdf(file = file.path(outDir, "ROC.genes.pdf"), width = 5.5, height = 5)
    plot(roc1, print.auc = FALSE, col = bioCol[k], legacy.axes = TRUE, main = "", lwd = 3)
    aucText <- c(aucText, paste0(x, ", AUC=", sprintf("%.3f", roc1$auc[1])))
  } else {
    plot(roc1, print.auc = FALSE, col = bioCol[k], legacy.axes = TRUE, main = "", lwd = 3, add = TRUE)
    aucText <- c(aucText, paste0(x, ", AUC=", sprintf("%.3f", roc1$auc[1])))
  }
}
legend("bottomright", aucText, lwd = 3, bty = "n", col = bioCol[1:(ncol(rt) - 1)])
dev.off()

# =====================================================================
# (1) PRIMARY: locked-coefficient extrapolation
# =====================================================================
expr <- rt[, 1:(ncol(rt) - 1), drop = FALSE]
locked.coef <- c(`(Intercept)` = 0.8263, ARG1 = -1.0028, GZMB = 0.3379, EPCAM = 0.1328)
lp <- locked.coef["(Intercept)"] +
  locked.coef["ARG1"] * expr$ARG1 +
  locked.coef["GZMB"] * expr$GZMB +
  locked.coef["EPCAM"] * expr$EPCAM
pred <- as.numeric(lp)          # linear predictor; AUC is invariant to the link transform

rocLocked <- roc(y, pred)
ciLocked <- ci.auc(rocLocked, method = "bootstrap")
ciVec <- as.numeric(ciLocked)
pdf(file = file.path(outDir, "ROC.model.locked.pdf"), width = 5.5, height = 5)
plot(rocLocked, print.auc = TRUE, col = "red", legacy.axes = TRUE,
     main = "Model (locked coefficients, external validation)", lwd = 3)
text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f", ciVec[1]), "-", sprintf("%.03f", ciVec[3])), col = "red")
dev.off()

bestLocked <- coords(rocLocked, x = "best", input = "threshold", best.method = "youden",
                     ret = c("threshold", "sensitivity", "specificity"))
cat("---- Locked-coefficient model (PRIMARY) ----\n")
print(ciLocked)
print(bestLocked)

# =====================================================================
# (2) SENSITIVITY: within-cohort refit (as in the original working script)
# =====================================================================
logit <- glm(y ~ ., family = binomial(link = "logit"), data = expr)
predRefit <- predict(logit, newdata = expr)    # FIXED: was predict(logit, newx=rt)

rocRefit <- roc(y, as.numeric(predRefit))
ciRefit <- ci.auc(rocRefit, method = "bootstrap")
ciVec2 <- as.numeric(ciRefit)
pdf(file = file.path(outDir, "ROC.model.refit.pdf"), width = 5.5, height = 5)
plot(rocRefit, print.auc = TRUE, col = "blue", legacy.axes = TRUE,
     main = "Model (within-cohort refit - SENSITIVITY ONLY)", lwd = 3)
text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f", ciVec2[1]), "-", sprintf("%.03f", ciVec2[3])), col = "blue")
dev.off()

bestRefit <- coords(rocRefit, x = "best", input = "threshold", best.method = "youden",
                    ret = c("threshold", "sensitivity", "specificity"))
cat("---- Within-cohort refit (SENSITIVITY) ----\n")
print(ciRefit)
print(bestRefit)
