# =====================================================================
# 06b - Venn intersection: DEGs x WGCNA turquoise module x NETs gene list
# Rebuilt from the analysis working file: "Version 1/9.Venn interGenes/NETs13.venn.R"
# Output: 7 intersection genes (ARG1, GZMB, EPCAM, IDO1, TNFRSF1B, CXCL1, CCL20)
#   = gene_lists/interGenes.txt (shipped)
# The NETs gene list was derived from a GeneCards "neutrophil extracellular
#   trap" search (the raw GeneCards export is not redistributed here for
#   licensing reasons; the derived symbol list is shipped in gene_lists/NETs.txt).
# Requires: results/04_deg/diff.txt (script 04);
#           gene_lists/module_turquoise.txt (= the turquoise module exported by script 06)
# =====================================================================

library(ggvenn)

outDir <- "results/06b_venn_interGenes"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

geneList <- list()

disRT <- read.table("results/04_deg/diff.txt", header = TRUE, sep = "\t", check.names = FALSE)
geneList[["DEGs"]] <- as.vector(disRT[, 1])

wgcnaRT <- read.table("gene_lists/module_turquoise.txt", header = FALSE, sep = "\t", check.names = FALSE)
geneList[["WGCNA"]] <- as.vector(wgcnaRT[, 1])

NETsRT <- read.table("gene_lists/NETs.txt", header = FALSE, sep = "\t", check.names = FALSE)
geneList[["NETs"]] <- as.vector(NETsRT[, 1])

pdf(file = file.path(outDir, "venn.pdf"), width = 6, height = 6)
ggvenn(geneList, show_percentage = TRUE, stroke_color = "white", stroke_size = 0.5,
       fill_color = c("#E41A1C", "#1E90FF", "#FF8C00"),
       set_name_color = c("#E41A1C", "#1E90FF", "#FF8C00"),
       set_name_size = 6, text_size = 4.5)
dev.off()

interGenes <- Reduce(intersect, geneList)
write.table(interGenes, file = file.path(outDir, "interGenes.txt"),
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
