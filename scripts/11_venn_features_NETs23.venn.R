# =====================================================================
# 11 - Venn intersection of the three feature-selection gene sets
# Rebuilt from the analysis working file: "Version 1/19.Venn interFeatureGenes/NETs23.venn.R"
# Output: interFeatureGenes.txt -> the 3-gene signature (ARG1, GZMB, EPCAM)
#   = gene_lists/interFeatureGenes.txt (shipped)
# Requires: results/08_lasso/LASSO.gene.txt, results/09_svm_rfe/SVM-RFE.gene.txt,
#           results/10_rf/RF.gene.txt
# =====================================================================

library(ggvenn)

outDir <- "results/11_venn_features"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

inputs <- c(file.path("results/08_lasso",  "LASSO.gene.txt"),
            file.path("results/09_svm_rfe", "SVM-RFE.gene.txt"),
            file.path("results/10_rf",     "RF.gene.txt"))

geneList <- list()
for (inputFile in inputs) {
  rt <- read.table(inputFile, header = FALSE, sep = "\t", check.names = FALSE, fill = TRUE)
  geneNames <- unlist(strsplit(as.vector(t(rt)), "\\s+"))
  geneNames <- geneNames[geneNames != "" & geneNames != "NA"]
  geneNames <- gsub("^ | $", "", geneNames)
  uniqGene <- unique(geneNames)
  setName <- gsub("\\.gene\\.txt$", "", basename(inputFile))
  geneList[[setName]] <- uniqGene
}

pdf(file = file.path(outDir, "venn.pdf"), width = 6, height = 6)
ggvenn(geneList, show_percentage = TRUE, stroke_color = "white", stroke_size = 0.5,
       fill_color = c("#E41A1C", "#1E90FF", "#FF8C00"),
       set_name_color = c("#E41A1C", "#1E90FF", "#FF8C00"),
       set_name_size = 6, text_size = 4.5)
dev.off()

interGenes <- Reduce(intersect, geneList)
write.table(interGenes, file = file.path(outDir, "interFeatureGenes.txt"),
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
