# =====================================================================
# 18a - Chromosomal location circle plot of the 7 intersection genes
# Rebuilt from the analysis working file: "Version 1/11.circle/NETs15.circlize.R"
# Requires: gene_lists/interGenes.txt, reference_data/geneREF.txt (shipped)
# =====================================================================

library(circlize)

outDir <- "results/18a_circlize"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

genepos <- read.table("reference_data/geneREF.txt", header = TRUE, sep = "\t", check.names = FALSE)
colnames(genepos) <- c("genename", "chr", "start", "end")
genepos <- genepos[, c("chr", "start", "end", "genename")]
row.names(genepos) <- genepos[, "genename"]

geneRT <- read.table("gene_lists/interGenes.txt", header = FALSE, sep = "\t", check.names = FALSE)
genepos <- genepos[as.vector(geneRT[, 1]), ]
bed0 <- genepos

pdf(file = file.path(outDir, "circlize.pdf"), width = 6, height = 6)
circos.clear()
circos.initializeWithIdeogram(species = "hg38", plotType = NULL)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr <- CELL_META$sector.index
  xlim <- CELL_META$xlim
  ylim <- CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(24))
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.6, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.15, bg.border = NA)
circos.genomicIdeogram(species = "hg38", track.height = mm_h(6))
circos.genomicLabels(bed0, labels.column = 4, side = "inside", cex = 0.8)
circos.clear()
dev.off()
