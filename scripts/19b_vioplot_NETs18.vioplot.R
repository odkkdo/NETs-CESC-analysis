# =====================================================================
# 19b - ssGSEA score heatmap and violin plots (Normal vs Cancer)
# Rebuilt from the analysis working file: "Version 1/14.vioplot/NETs18.vioplot.R"
# MINOR FIX: the zero-variance guard wrote the nudge to the wrong row
#   (rt[(conNum+1),i] for the tumor block); corrected to a random index of the
#   respective group. No effect on results unless a group SD is exactly 0.
# Requires: results/07_ssGSEA/ssGSEA.result.txt
# =====================================================================

library(vioplot)
library(pheatmap)

outDir <- "results/19b_vioplot"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

rt <- read.table("results/07_ssGSEA/ssGSEA.result.txt", header = TRUE, sep = "\t",
                 check.names = FALSE, row.names = 1)

con   <- grepl("_Normal", colnames(rt), ignore.case = TRUE)
treat <- grepl("_Cancer", colnames(rt), ignore.case = TRUE)
conData <- rt[, con]
treatData <- rt[, treat]
conNum <- ncol(conData)
treatNum <- ncol(treatData)
rt <- t(cbind(conData, treatData))

# ---- heatmap ----
Type <- c(rep("Normal", conNum), rep("Cancer", treatNum))
names(Type) <- row.names(rt)
Type <- as.data.frame(Type)
pdf(file = file.path(outDir, "heatmap.pdf"), width = 8, height = 5)
pheatmap(t(rt),
         annotation_col = Type,
         color = colorRampPalette(c(rep("#4DBBD5", 3), "white", rep("#E64B35", 3)))(50),
         cluster_cols = FALSE,
         show_colnames = FALSE,
         scale = "row",
         fontsize = 6,
         fontsize_row = 6,
         fontsize_col = 6)
dev.off()

# ---- violin plots with Wilcoxon p values ----
pdf(file = file.path(outDir, "vioplot.pdf"), width = 12, height = 7.5)
par(las = 1, mar = c(10, 6, 3, 3))
x <- c(1:ncol(rt)); y <- c(1:ncol(rt))
plot(x, y, xlim = c(0, (3 * ncol(rt) - 3)), ylim = c(min(rt), max(rt) + 0.05),
     main = "", xlab = "", ylab = "Fraction", pch = 21, col = "white", xaxt = "n")
for (i in 1:ncol(rt)) {
  # avoid zero SD causing vioplot to fail (nudge a random sample of that group)
  if (sd(rt[1:conNum, i]) == 0) { rt[sample(1:conNum, 1), i] <- 0.001 }
  if (sd(rt[(conNum + 1):(conNum + treatNum), i]) == 0) {
    rt[sample((conNum + 1):(conNum + treatNum), 1), i] <- 0.001
  }

  conData <- rt[1:conNum, i]
  treatData <- rt[(conNum + 1):(conNum + treatNum), i]

  vioplot(conData, at = 3 * (i - 1), col = "#4DBBD5", border = "black", lty = 1, add = TRUE)
  vioplot(treatData, at = 3 * (i - 1) + 1, col = "#E64B35", border = "black", lty = 1, add = TRUE)

  wilcoxTest <- wilcox.test(conData, treatData)
  p <- wilcoxTest$p.value

  mx <- max(c(conData, treatData))
  lines(c(3 * (i - 1) + 0.1, 3 * (i - 1) + 0.9), c(mx, mx), lwd = 1.2)
  text(x = 3 * (i - 1) + 0.5, y = mx + 0.03,
       labels = ifelse(p < 0.001, "p<0.001", paste0("p=", sprintf("%.03f", p))),
       cex = 0.8, col = "black")
}

legend("topright", legend = c("Normal", "Cancer"), fill = c("#4DBBD5", "#E64B35"),
       border = "black", bty = "n")
text(seq(1, (3 * ncol(rt) - 2), 3), -0.075, xpd = NA, labels = colnames(rt),
     cex = 1, srt = 45, pos = 2)
dev.off()
