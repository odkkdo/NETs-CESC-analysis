# =====================================================================
# 13 - Nomogram, calibration curve and decision-curve analysis
# Rebuilt from the analysis working file: "Version 1/22.Nomo/Nomo.R"
# (paths made repo-relative; no methodological changes)
# Requires: data/merge.normalize.txt, gene_lists/interFeatureGenes.txt
# =====================================================================

library(rms)
library(rmda)

outDir <- "results/13_nomogram"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

data <- read.table("data/merge.normalize.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
geneRT <- read.table("gene_lists/interFeatureGenes.txt", header = FALSE, sep = "\t", check.names = FALSE)
data <- data[as.vector(geneRT[, 1]), ]

data <- t(data)
group <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(data))
rt <- cbind(as.data.frame(data), Type = group)
rt$Type <- factor(rt$Type, levels = c("Normal", "Cancer"))

ddist <- datadist(rt)
options(datadist = "ddist")

lrmModel <- lrm(Type ~ ., data = rt, x = TRUE, y = TRUE)
nomo <- nomogram(lrmModel, fun = plogis,
                 fun.at = c(0.001, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99),
                 lp = FALSE, funlabel = "Risk of Disease")
pdf(file = file.path(outDir, "Nomo.pdf"), width = 8, height = 6)
plot(nomo)
dev.off()

cali <- calibrate(lrmModel, method = "boot", B = 1000)
pdf(file = file.path(outDir, "Calibration.pdf"), width = 5.5, height = 5.5)
plot(cali, xlab = "Predicted probability", ylab = "Actual probability", sub = FALSE)
dev.off()

rt$Type <- ifelse(rt$Type == "Normal", 0, 1)
dc <- decision_curve(Type ~ ARG1 + GZMB + EPCAM, data = rt,
                     family = binomial(link = "logit"),
                     thresholds = seq(0, 1, by = 0.01),
                     confidence.intervals = 0.95)
pdf(file = file.path(outDir, "DCA.pdf"), width = 6, height = 6)
plot_decision_curve(dc, curve.names = "Model", xlab = "Threshold probability",
                    cost.benefit.axis = TRUE, col = "red",
                    confidence.intervals = FALSE, standardize = FALSE)
dev.off()
