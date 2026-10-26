# =====================================================================
# 03 - PCA before / after batch correction
# Rebuilt from the analysis working file: "Version 1/4.PCA/geo07.PCA.R"
# (default argument `null` corrected to NULL; paths made repo-relative)
# Run from the repository root.
# =====================================================================

library(ggplot2)
library(ggpubr)

outDir <- "results/03_PCA"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

bioPCA <- function(inputFile = NULL, outFile = NULL, titleName = NULL) {
  rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  data <- t(rt)
  Project <- gsub("(.*?)\\_.*", "\\1", rownames(data))   # GEO study id prefix

  data.pca  <- prcomp(data)
  pcaPredict <- predict(data.pca)
  PCA <- data.frame(PC1 = pcaPredict[, 1], PC2 = pcaPredict[, 2], Type = Project)

  pdf(file = outFile, width = 5.5, height = 4.25)
  p1 <- ggscatter(data = PCA, x = "PC1", y = "PC2", color = "Type", shape = "Type",
                  ellipse = TRUE, ellipse.type = "norm", ellipse.border.remove = FALSE,
                  ellipse.alpha = 0.1, size = 2, main = titleName, legend = "right") +
        theme(plot.margin = unit(rep(1.5, 4), "lines"),
              plot.title = element_text(hjust = 0.5))
  print(p1)
  dev.off()
}

bioPCA(inputFile = "data/merge.preNorm.txt",   # produced by script 02 (not shipped)
       outFile   = file.path(outDir, "PCA.preNorm.pdf"),
       titleName = "Before batch correction")
bioPCA(inputFile = "data/merge.normalize.txt", # shipped in data/
       outFile   = file.path(outDir, "PCA.normalize.pdf"),
       titleName = "After batch correction")
