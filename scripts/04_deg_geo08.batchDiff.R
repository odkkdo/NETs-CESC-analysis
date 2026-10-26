# =====================================================================
# 04 - Differential expression analysis (limma), Cancer vs Normal
# Rebuilt from the analysis working file: "Version 1/5.差异基因/geo08.batchDiff.R"
#
# METHOD FIX (documented in the response letter, applied to the revision):
#   The original design matrix `model.matrix(~0 + factor(Type))` ordered the
#   factor levels alphabetically (Cancer, Normal) but then renamed the columns
#   c("Normal","Cancer"). The contrast "Cancer - Normal" therefore actually
#   computed Normal - Cancer: every logFC written to diff.txt / all.txt (and
#   the Up/Down colours of the volcano plot) had the wrong sign.
#   Fix: set the factor levels explicitly (Normal, Cancer) and drop the manual
#   column renaming. Biological conclusions are unchanged (ARG1 down-regulated
#   in tumors, the other six genes up-regulated; manuscript Table 2 already
#   used the correct directions).
# Cosmetic repairs: pheapheatmap( -> pheatmap(,  def() -> dev.off(), and the
#   broken assignment line before the volcano plot restored.
# =====================================================================

library(limma)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggrepel)

logFCfilter      <- 1      # |log2FC| > 1 (2-fold)
adj.P.Val.Filter <- 0.05   # Benjamini-Hochberg adjusted P
inputFile        <- "data/merge.normalize.txt"
outDir <- "results/04_deg"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

# ---- read expression matrix ----
rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)

# ---- sample annotations (sample names: <GEO study>_<GSM id>_<Type>) ----
Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data <- data[, order(Type)]
Project <- gsub("(.+)\\_(.+)\\_(.+)", "\\1", colnames(data))
Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
colnames(data) <- gsub("(.+)\\_(.+)\\_(.+)", "\\2", colnames(data))

# ---- limma model (FIXED: explicit factor levels) ----
design <- model.matrix(~ 0 + factor(Type, levels = c("Normal", "Cancer")))
colnames(design) <- c("Normal", "Cancer")
fit  <- lmFit(data, design)
cont.matrix <- makeContrasts(Cancer - Normal, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# ---- all genes ----
allDiff <- topTable(fit2, adjust = "fdr", number = 200000)
allDiffOut <- rbind(id = colnames(allDiff), allDiff)
write.table(allDiffOut, file = file.path(outDir, "all.txt"), sep = "\t", quote = FALSE, col.names = FALSE)

# ---- significant DEGs: |logFC| > 1 and BH adj.P < 0.05 (123 DEGs in the study) ----
diffSig <- allDiff[with(allDiff, (abs(logFC) > logFCfilter & adj.P.Val < adj.P.Val.Filter)), ]
diffSigOut <- rbind(id = colnames(diffSig), diffSig)
write.table(diffSigOut, file = file.path(outDir, "diff.txt"), sep = "\t", quote = FALSE, col.names = FALSE)

diffGeneExp <- data[row.names(diffSig), ]
diffGeneExpOut <- rbind(id = paste0(colnames(diffGeneExp), "_", Type), diffGeneExp)
write.table(diffGeneExpOut, file = file.path(outDir, "diffGeneExp.txt"), sep = "\t", quote = FALSE, col.names = FALSE)

# ---- heatmap of top DEGs ----
geneNum <- 50
diffUp   <- diffSig[diffSig$logFC > 0, ]
diffDown <- diffSig[diffSig$logFC < 0, ]
geneUp   <- row.names(diffUp)
geneDown <- row.names(diffDown)
if (nrow(diffUp)   > geneNum) { geneUp   <- row.names(diffUp)[1:geneNum] }
if (nrow(diffDown) > geneNum) { geneDown <- row.names(diffDown)[1:geneNum] }
hmExp <- data[c(geneUp, geneDown), ]
names(Type) <- colnames(data)
Type <- as.data.frame(Type)
Type <- cbind(Project, Type)
pdf(file = file.path(outDir, "heatmap.pdf"), width = 10, height = 7)
pheatmap(hmExp,                              # FIXED: was "pheapheatmap("
         annotation_col = Type,
         color = colorRampPalette(c("blue2", "white", "red2"))(50),
         cluster_cols = FALSE,
         show_colnames = FALSE,
         scale = "row",
         fontsize = 8,
         fontsize_row = 5.5,
         fontsize_col = 8)
dev.off()                                    # FIXED: was "def()"

# ---- volcano plot (Fig. 2A; Up/Down now follows the corrected contrast) ----
rt <- read.table(file.path(outDir, "all.txt"), header = TRUE, sep = "\t", check.names = FALSE)
Sig <- ifelse((rt$adj.P.Val < adj.P.Val.Filter) & (abs(rt$logFC) > logFCfilter),
              ifelse(rt$logFC > logFCfilter, "Up", "Down"), "Not")
rt <- mutate(rt, Sig = Sig)

topUp <- rt %>% filter(Sig == "Up")   %>% arrange(adj.P.Val) %>% head(5)
topDown <- rt %>% filter(Sig == "Down") %>% arrange(adj.P.Val) %>% head(5)
topGene <- rbind(topUp, topDown)

p <- ggplot(rt, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Sig), size = 1.6, alpha = 0.8) +
  scale_color_manual(values = c("Down" = "#3C8DBC", "Not" = "grey80", "Up" = "#E64B35")) +
  geom_vline(xintercept = c(-logFCfilter, logFCfilter), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(adj.P.Val.Filter), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = topGene, aes(label = topGene$id), size = 4,
                  box.padding = 0.4, point.padding = 0.3,
                  segment.color = "grey40", max.overlaps = 100) +
  labs(x = "log2 Fold Change", y = "-log10(adj.P.Val)", title = "Volcano Plot", color = "") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        panel.grid = element_blank())
pdf(file = file.path(outDir, "vol.pdf"), width = 5.5, height = 4.5)
print(p)
dev.off()
