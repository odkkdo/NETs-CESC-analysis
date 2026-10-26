# =====================================================================
# 12 - Training-cohort ROC: single genes + logistic model
# Rebuilt from the analysis working file: "Version 1/21.ROC/NETs25.ROC.R"
# FIXES vs the working file:
#   1) y=ifelse(y=="Control",0,1) -> y=="Normal" (the cohort labels are
#      "Normal"/"Cancer", so the original line labelled every sample as case);
#   2) pd0000f( -> pdf(;
#   3) predict(logit, newx=rt) -> predict(logit, newdata=rt): `newx` is not a
#      predict.glm argument. The model ROC is an in-sample (training-cohort)
#      estimate; the reported training AUC = 0.864 (95% CI 0.770-0.934,
#      bootstrap) refers to this model. External validation is handled by
#      script 14 (locked-coefficient extrapolation).
# Requires: data/merge.normalize.txt, gene_lists/interFeatureGenes.txt
# =====================================================================

library(glmnet)
library(pROC)

outDir <- "results/12_ROC"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

rt <- read.table("data/merge.normalize.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

y <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(rt))
y <- ifelse(y == "Normal", 0, 1)          # FIXED: was y=="Control"

geneRT <- read.table("gene_lists/interFeatureGenes.txt", header = FALSE, sep = "\t", check.names = FALSE)

bioCol <- c("#F05C3BFF", "#5C88DAFF", "#5CB85CFF", "#EEA236FF", "#9632B8FF", "#17BECFFF", "#BCBD22FF")
bioCol <- bioCol[1:nrow(geneRT)]

# ---- per-gene ROC curves ----
aucText <- c()
k <- 0
for (x in as.vector(geneRT[, 1])) {
  k <- k + 1
  roc1 <- roc(y, as.numeric(rt[x, ]))
  if (k == 1) {
    pdf(file = file.path(outDir, "ROC.genes.pdf"), width = 5.5, height = 5)
    plot(roc1, print.auc = FALSE, col = bioCol[k], legacy.axes = TRUE, main = "", lwd = 3)
    aucText <- c(aucText, paste0(x, ", AUC=", sprintf("%.3f", roc1$auc[1])))
  } else {
    plot(roc1, print.auc = FALSE, col = bioCol[k], legacy.axes = TRUE, main = "", add = TRUE, lwd = 3)
    aucText <- c(aucText, paste0(x, ", AUC=", sprintf("%.3f", roc1$auc[1])))
  }
}
legend("bottomright", aucText, lwd = 3, bty = "n", cex = 1, col = bioCol[1:length(aucText)])
dev.off()

# ---- logistic model on the 3 signature genes (in-sample ROC) ----
rt2 <- as.data.frame(t(rt[as.vector(geneRT[, 1]), ]))
logit <- glm(y ~ ., family = binomial(link = "logit"), data = rt2)
pred <- predict(logit, newdata = rt2)     # FIXED: was newx=rt

roc1 <- roc(y, as.numeric(pred))
ci1 <- ci.auc(roc1, method = "bootstrap")
ciVec <- as.numeric(ci1)
pdf(file = file.path(outDir, "ROC.model.pdf"), width = 5.5, height = 5)
plot(roc1, print.auc = TRUE, col = "red", legacy.axes = TRUE, main = "Model", lwd = 3)
text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f", ciVec[1]), "-", sprintf("%.03f", ciVec[3])), col = "red")
dev.off()
