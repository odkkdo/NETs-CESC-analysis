# =====================================================================
# rebuild_fig2A_and_5J.R
# Regenerate Figure 2A (volcano plot, corrected DEG direction) and
# Figure 5J (external-validation ROC, locked coefficients as primary)
# with the two methodological fixes documented in the response letter.
#
# Outputs (PNG 300 dpi + PDF) go to OUT_DIR, default "E:/part-time/WWW/New figure":
#   fig2A_volcano.png / .pdf   -> replace panel A of manuscript Figure 2 (1.tif)
#   fig5J_ROC.png / .pdf       -> replace panel J of manuscript Figure 5 (4.tif)
#
# Requirements: R >= 4.2 with packages limma, ggplot2, ggrepel, pROC.
# Run:  Rscript rebuild_fig2A_and_5J.R
# =====================================================================

REPO    <- "C:/Users/Administrator/AppData/Roaming/AutoClaw-official/accounts/b1d764ea8ea052d069a7b085230b63a2db22616de0ed058e8a24a891f8c9040b/projects/Default/NETs-CESC-analysis"
EXT_MAT <- "F:/WWW/3-4trik/NETs-CESC/Version 1/23.testDiff/GTExTCGA.merge.txt"  # 162 MB external matrix
OUT_DIR <- "E:/part-time/WWW/New figure"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

library(limma)
library(ggplot2)
library(ggrepel)
library(pROC)

## ------------------------------------------------------------------ ##
## Figure 2A - volcano plot from the CORRECTED DEG analysis            ##
## (design matrix factor levels fixed: Normal, Cancer)                 ##
## ------------------------------------------------------------------ ##
rt <- read.table(file.path(REPO, "data/merge.normalize.txt"), header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)

Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data <- data[, order(Type)]
Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
colnames(data) <- gsub("(.+)\\_(.+)\\_(.+)", "\\2", colnames(data))

# FIXED design: explicit factor levels (original working file computed Normal - Cancer)
design <- model.matrix(~ 0 + factor(Type, levels = c("Normal", "Cancer")))
colnames(design) <- c("Normal", "Cancer")
fit  <- lmFit(data, design)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(Cancer - Normal, levels = design)))

allDiff <- topTable(fit2, adjust = "fdr", number = 200000)
logFCfilter <- 1; adj.P.Val.Filter <- 0.05
allDiff$Sig <- ifelse(allDiff$adj.P.Val < adj.P.Val.Filter & abs(allDiff$logFC) > logFCfilter,
                      ifelse(allDiff$logFC > 0, "Up", "Down"), "Not")
cat("DEGs (|logFC|>1, BH<0.05):", sum(allDiff$Sig != "Not"),
    "| Up:", sum(allDiff$Sig == "Up"), "| Down:", sum(allDiff$Sig == "Down"), "\n")
cat("ARG1 logFC (expect negative):", allDiff["ARG1", "logFC"], "\n")
cat("GZMB logFC (expect positive):", allDiff["GZMB", "logFC"], "\n")

topGene <- rbind(head(allDiff[allDiff$Sig == "Up" & !rownames(allDiff) %in% c("ARG1"), ], 5),
                 head(allDiff[allDiff$Sig == "Down", ], 5))
topGene$gene <- rownames(topGene)

p2a <- ggplot(allDiff, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Sig), size = 1.6, alpha = 0.8) +
  scale_color_manual(values = c("Down" = "#3C8DBC", "Not" = "grey80", "Up" = "#E64B35"),
                     breaks = c("Up", "Not", "Down")) +
  geom_vline(xintercept = c(-logFCfilter, logFCfilter), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(adj.P.Val.Filter), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = topGene, aes(label = gene), size = 3.2,
                  box.padding = 0.4, point.padding = 0.3,
                  segment.color = "grey40", max.overlaps = 100) +
  labs(x = "log2 Fold Change (tumor vs. normal)", y = "-log10 (adjusted P)",
       title = NULL, color = "") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        legend.position = "right")

ggsave(file.path(OUT_DIR, "fig2A_volcano.png"), p2a, width = 5.5, height = 4.5, dpi = 300)
ggsave(file.path(OUT_DIR, "fig2A_volcano.pdf"), p2a, width = 5.5, height = 4.5)

## ------------------------------------------------------------------ ##
## Figure 5J - external validation ROC (n = 404)                       ##
## PRIMARY: locked training coefficients applied to external matrix    ##
## SENSITIVITY: within-cohort refit (as in the original report)        ##
## ------------------------------------------------------------------ ##
ext <- read.table(EXT_MAT, header = TRUE, sep = "\t", check.names = FALSE)
ext <- as.matrix(ext)
rownames(ext) <- ext[, 1]
ex <- ext[, 2:ncol(ext)]
dimnames <- list(rownames(ex), colnames(ex))
edata <- matrix(as.numeric(as.matrix(ex)), nrow = nrow(ex), dimnames = dimnames)
edata <- avereps(edata)

genes3 <- c("ARG1", "GZMB", "EPCAM")
stopifnot(all(genes3 %in% rownames(edata)))
em <- t(edata[genes3, ])
TypeE <- gsub("(.*)\\_(.*)", "\\2", rownames(em))
y <- ifelse(TypeE == "Normal", 0, 1)
cat("External cohort: n =", length(y), "| tumors:", sum(y == 1), "| normals:", sum(y == 0), "\n")

# PRIMARY: locked coefficients from the 78-sample training cohort
locked <- c(`(Intercept)` = 0.8263, ARG1 = -1.0028, GZMB = 0.3379, EPCAM = 0.1328)
lp <- locked["(Intercept)"] + locked["ARG1"] * em[, "ARG1"] +
  locked["GZMB"] * em[, "GZMB"] + locked["EPCAM"] * em[, "EPCAM"]
rocLocked <- roc(y, as.numeric(lp), levels = c(0, 1), direction = "<")
ciL <- as.numeric(ci.auc(rocLocked, method = "bootstrap"))
bL <- coords(rocLocked, x = "best", input = "threshold", best.method = "youden",
             ret = c("threshold", "sensitivity", "specificity"))
cat(sprintf("LOCKED  AUC=%.3f (95%% CI %.3f-%.3f) | Youden thr=%.3f Sen=%.3f Spec=%.3f\n",
            rocLocked$auc[1], ciL[1], ciL[3], bL$threshold, bL$sensitivity, bL$specificity))

# SENSITIVITY: within-cohort refit
dfE <- as.data.frame(em)
refit <- glm(y ~ ., family = binomial(link = "logit"), data = dfE)
predR <- as.numeric(predict(refit, newdata = dfE, type = "link"))
rocRefit <- roc(y, predR, levels = c(0, 1), direction = "<")
ciR <- as.numeric(ci.auc(rocRefit, method = "bootstrap"))
bR <- coords(rocRefit, x = "best", input = "threshold", best.method = "youden",
             ret = c("threshold", "sensitivity", "specificity"))
cat(sprintf("REFIT   AUC=%.3f (95%% CI %.3f-%.3f) | Youden thr=%.3f Sen=%.3f Spec=%.3f\n",
            rocRefit$auc[1], ciR[1], ciR[3], bR$threshold, bR$sensitivity, bR$specificity))

rocData <- rbind(
  data.frame(fpr = 1 - as.numeric(rocLocked$specificities),
             tpr = as.numeric(rocLocked$sensitivities),
             Model = sprintf("Locked coefficients (AUC=%.3f, 95%% CI %.3f-%.3f)",
                             rocLocked$auc[1], ciL[1], ciL[3])),
  data.frame(fpr = 1 - as.numeric(rocRefit$specificities),
             tpr = as.numeric(rocRefit$sensitivities),
             Model = sprintf("Within-cohort refit (AUC=%.3f, 95%% CI %.3f-%.3f)",
                             rocRefit$auc[1], ciR[1], ciR[3])))

rocData$Model <- factor(rocData$Model, levels = c(
  sprintf("Locked coefficients (AUC=%.3f, 95%% CI %.3f-%.3f)", rocLocked$auc[1], ciL[1], ciL[3]),
  sprintf("Within-cohort refit (AUC=%.3f, 95%% CI %.3f-%.3f)", rocRefit$auc[1], ciR[1], ciR[3])))

p5j <- ggplot(rocData, aes(x = fpr, y = tpr, color = Model)) +
  geom_path(lwd = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = setNames(c("#E64B35", "#1F4E79"), levels(rocData$Model))) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = NULL, color = NULL) +
  coord_equal() +
  theme_bw(base_size = 13) +
  theme(legend.position = c(0.62, 0.22),
        legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
        legend.text = element_text(size = 9),
        panel.grid = element_blank())

ggsave(file.path(OUT_DIR, "fig5J_ROC.png"), p5j, width = 5.5, height = 5.0, dpi = 300)
ggsave(file.path(OUT_DIR, "fig5J_ROC.pdf"), p5j, width = 5.5, height = 5.0)

writeLines(c(
  "Figure regeneration summary",
  "",
  sprintf("Fig 2A: DEGs Up=%d Down=%d; ARG1 logFC=%.3f (down), GZMB logFC=%.3f (up)",
          sum(allDiff$Sig == "Up"), sum(allDiff$Sig == "Down"),
          allDiff["ARG1", "logFC"], allDiff["GZMB", "logFC"]),
  sprintf("Fig 5J locked : AUC=%.3f CI %.3f-%.3f | Youden %.3f -> Sen %.1f%% Spec %.1f%%",
          rocLocked$auc[1], ciL[1], ciL[3], bL$threshold, 100 * bL$sensitivity, 100 * bL$specificity),
  sprintf("Fig 5J refit  : AUC=%.3f CI %.3f-%.3f | Youden %.3f -> Sen %.1f%% Spec %.1f%%",
          rocRefit$auc[1], ciR[1], ciR[3], bR$threshold, 100 * bR$sensitivity, 100 * bR$specificity)),
  file.path(OUT_DIR, "figure_summary.txt"))

cat("Done. Outputs in:", OUT_DIR, "\n")
