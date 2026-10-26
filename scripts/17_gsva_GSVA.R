# =====================================================================
# 17 - GSVA metabolic pathway analysis per molecular subtype
# Rebuilt from the analysis working file: "Version 1/28.GSVA/GSVA.R"
# Parameters: KEGG gene sets (MSigDB c2.cp.kegg symbols), gsvaParam
#   minSize = 5 / maxSize = 500; limma per-subtype vs mean-of-others
#   contrasts. In the study: C2 cytosolic DNA sensing (BH adj.P = 0.043) and
#   C3 gap junction (BH adj.P = 0.021) survived correction; the C1 top-3
#   pathways (pyruvate metabolism / ubiquitin-mediated proteolysis /
#   melanogenesis) were nominal trends only.
# Requires: data/merge.normalize.txt, results/15_cluster/cluster.txt and a
#   KEGG gene-set file reference_data/c2.cp.kegg.Hs.symbols.gmt (download
#   from MSigDB; not redistributed here).
# =====================================================================

library(reshape2)
library(ggpubr)
library(limma)
library(dplyr)
library(tidyr)
library(GSEABase)
library(GSVA)
library(ggplot2)

expFile <- "data/merge.normalize.txt"
clusterFile <- "results/15_cluster/cluster.txt"
gmtFile <- "reference_data/c2.cp.kegg.Hs.symbols.gmt"
outDir <- "results/17_GSVA"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. expression matrix, tumors only ----
rt <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)

group <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data <- data[, group == "Cancer", drop = FALSE]

# ---- 2. GSVA ----
geneSets <- getGmt(gmtFile, geneIdType = SymbolIdentifier())
gsvaPar <- gsvaParam(data, geneSets, minSize = 5, maxSize = 500)
gsvaScore <- gsva(gsvaPar)

normalize <- function(x) { (x - min(x)) / (max(x) - min(x)) }
gsvaScore <- normalize(gsvaScore)

# ---- 3. cluster assignments ----
cluster <- read.table(clusterFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
groups <- unique(cluster$Cluster)   # "C1", "C2", "C3"

sampleList <- lapply(groups, function(g) rownames(cluster)[cluster$Cluster == g])
names(sampleList) <- groups

dataList <- lapply(sampleList, function(s) gsvaScore[, s, drop = FALSE])
data_combined <- do.call(cbind, dataList)

Type <- factor(rep(groups, times = sapply(sampleList, length)))

# ---- 4. limma per-subtype vs mean-of-others ----
design <- model.matrix(~ 0 + Type)
colnames(design) <- levels(Type)
fit <- lmFit(data_combined, design)

contrast.matrix <- makeContrasts(
  C1vsOthers = C1 - (C2 + C3) / 2,
  C2vsOthers = C2 - (C1 + C3) / 2,
  C3vsOthers = C3 - (C1 + C2) / 2,
  levels = design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

top3_list <- list()
for (comp in colnames(contrast.matrix)) {
  tt <- topTable(fit2, coef = comp, number = Inf, sort.by = "P")
  tt$Pathway <- rownames(tt)
  tt$Comparison <- comp
  tt_up <- tt[tt$logFC > 0, ]
  if (nrow(tt_up) == 0) {
    warning(paste("No logFC > 0 pathway in", comp, "- skipped"))
    next
  }
  top3 <- tt_up[order(tt_up$P.Value), ][1:min(3, nrow(tt_up)), ]
  top3_list[[comp]] <- top3
}
top3_all <- do.call(rbind, top3_list)
top3_all$Cluster <- gsub("vsOthers", "", top3_all$Comparison)

print(top3_all[, c("Pathway", "Cluster", "logFC", "P.Value", "adj.P.Val")])
if (nrow(top3_all) == 0) { stop("No qualifying pathways found.") }

# ---- 5. plotting data ----
pathway_list <- unique(top3_all$Pathway)
plot_df <- expand.grid(Pathway = pathway_list, Cluster = levels(Type), stringsAsFactors = FALSE)

for (i in 1:nrow(plot_df)) {
  p <- plot_df$Pathway[i]
  cl <- plot_df$Cluster[i]
  idx <- which(Type == cl)
  expr_vals <- as.numeric(data_combined[p, idx])
  plot_df$Mean[i] <- mean(expr_vals)
  plot_df$SD[i] <- sd(expr_vals)
}

belongs <- top3_all[, c("Pathway", "Cluster")]
names(belongs)[2] <- "BelongsTo"
plot_df <- left_join(plot_df, belongs, by = "Pathway")
plot_df$Pathway_label <- paste0(plot_df$BelongsTo, ": ", plot_df$Pathway)

belong_means <- plot_df %>%
  filter(Cluster == BelongsTo) %>%
  select(Pathway, BelongMean = Mean)
plot_df <- plot_df %>%
  left_join(belong_means, by = "Pathway") %>%
  arrange(factor(BelongsTo, levels = c("C1", "C2", "C3")), desc(BelongMean)) %>%
  mutate(Pathway_label = factor(Pathway_label, levels = unique(Pathway_label)))

# ---- 6. horizontal bar chart ----
cluster_colors <- c("#0066FF", "#FF0000", "#FF9900")
names(cluster_colors) <- levels(Type)

p <- ggplot(plot_df, aes(x = Mean, y = Pathway_label, fill = Cluster)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(xmin = Mean - SD, xmax = Mean + SD),
                position = position_dodge(0.8), width = 0.2, linewidth = 0.4) +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Top 3 Upregulated Pathways per Cluster",
       x = "GSVA Score (Mean \u00b1 SD)", y = NULL, fill = "Cluster") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 11, color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 12),
        plot.margin = margin(10, 10, 10, 10)) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15)))

pdf(file.path(outDir, "top3_upregulated_pathways_per_cluster_allbars.pdf"),
    width = 14, height = max(6, length(pathway_list) * 0.8))
print(p)
dev.off()

# ---- 7. result table ----
write.table(top3_all, file = file.path(outDir, "top3_upregulated_pathways_per_cluster.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
