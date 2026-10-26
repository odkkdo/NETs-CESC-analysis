# NETs-CESC-analysis

R scripts for the study:

> **A Three-Gene Neutrophil Extracellular Traps-Related Signature Enables Diagnosis and Immune-Metabolic Molecular Subtyping of Cervical Cancer**
> (manuscript under review; three-gene signature **ARG1 / GZMB / EPCAM**)

The pipeline integrates two GEO microarray/NanoString cohorts into a 78-sample
discovery cohort (38 normal + 40 cervical cancer), derives a NETs-related gene
signature, builds and validates a diagnostic model, and performs
immune-metabolic molecular subtyping of the tumors.

---

## 1. Pipeline overview

Scripts are numbered in execution order and are run from the repository root
(`setwd("<repo-root>")` before sourcing). Outputs are written to `results/`.

| Script | Step | Key inputs (not produced in-repo unless noted) | Key outputs |
|---|---|---|---|
| `scripts/01_normalize_GEOnormolize.R` | Per-dataset normalization (limma) | raw `geneMatrix.txt`, sample lists per GEO dataset | `<GSE>.normalize.txt` |
| `scripts/02_merge_combat_geo06.sva.R` | Merge datasets + ComBat batch correction | outputs of 01 | `merge.normalize.txt` (shipped in `data/`) |
| `scripts/03_pca_geo07.PCA.R` | PCA before/after batch correction | 02 outputs | PCA plots |
| `scripts/04_deg_geo08.batchDiff.R` | **DEG analysis (limma), FIXED** | `data/merge.normalize.txt` | 123 DEGs (`diff.txt`), volcano, heatmap |
| `scripts/05a/05b` | GO / KEGG enrichment of all DEGs | 04 output | enrichment tables & plots |
| `scripts/06_wgcna_NETs11.WGCNA.R` | WGCNA (beta = 8) | `data/merge.normalize.txt` | modules; turquoise = 148 genes |
| `scripts/06b_venn_interGenes_NETs13.venn.R` | DEGs ∩ turquoise ∩ NETs list | 04, 06, `gene_lists/NETs.txt` | 7 genes (`interGenes.txt`) |
| `scripts/07_ssgsea_NETs17.ssGSEA.R` | ssGSEA, 28-signature panel (22 scoreable) | `reference_data/immune.gmt` | ssGSEA scores |
| `scripts/08_lasso_NETs20.lasso.R` | LASSO feature selection (seed added) | 7 genes | 3 genes |
| `scripts/09_svm_rfe_NETs21.SVM.R` | SVM-RFE (caret, seed added) | 7 genes | ranked features |
| `scripts/10_rf_NETs22.randomForest.R` | Random forest importance (seed 12345) | 7 genes | top features |
| `scripts/11_venn_features_NETs23.venn.R` | LASSO ∩ SVM-RFE ∩ RF | 08/09/10 outputs | **ARG1, GZMB, EPCAM** |
| `scripts/12_roc_NETs25.ROC.R` | Training-cohort ROC (**FIXED**) | 3 genes | training AUC = 0.864 (95% CI 0.770–0.934) |
| `scripts/13_nomogram_Nomo.R` | Nomogram, calibration, DCA | 3 genes | nomogram/Calibration/DCA |
| `scripts/14_external_validation_testDiff.R` | **External validation, FIXED** | `data/GTExTCGA.merge.txt` (see below) | locked model AUC = 0.968 (0.940–0.989); refit AUC = 0.994 (sensitivity) |
| `scripts/15_clustering_Cluster.R` | ConsensusClusterPlus, tumors only (K = 3) | 3 genes | C1 = 20, C2 = 12, C3 = 8 |
| `scripts/16_subtype_immune_clusterlmm.R` | Immune scores per subtype | 15, 07 | immune.diff plot |
| `scripts/17_gsva_GSVA.R` | GSVA KEGG pathways per subtype | 15, MSigDB c2.cp.kegg (user download) | top pathways |
| `scripts/18a/18b` | Chromosome circle plot / Manhattan plot | `reference_data/geneREF.txt` | plots |
| `scripts/19a/19b/19c` | Expression boxplots / violin plots / gene×immune Spearman correlations (tumors only) | 07, 04 outputs | plots + correlation table |
| `scripts/20_drug_enrichment_NETs34.drugEnrich.R` | Drug enrichment (DSigDB) | 3 genes, DSigDB table (user download) | 14 candidate drugs |
| `scripts/21_sessionInfo_template.R` | Record R session | — | `sessionInfo.txt` |

## 2. Environment

* R ≥ 4.2.0 (the original analysis used R 4.2.0; run `scripts/21_sessionInfo_template.R` to record yours).
* Packages: `limma, sva, WGCNA, GSVA, GSEABase, clusterProfiler, org.Hs.eg.db, enrichplot, ComplexHeatmap, ConsensusClusterPlus, glmnet, randomForest, caret, e1071, pROC, rms, rmda, circlize, ggplot2, pheatmap, ggpubr, reshape2, CMplot, ggvenn, corrplot, vioplot, dplyr, tidyr, tidyverse, tidytable, RColorBrewer, ggrepel`.

## 3. Data sources

| Cohort | Source | Platform | Retained samples |
|---|---|---|---|
| GSE63514 | GEO (GSE63514) | GPL570 (Affymetrix HG-U133 Plus 2.0) | 24 normal + 28 cancer (full series: 24 normal + 76 CIN [14 low-, 22 intermediate-, 40 high-grade] + 28 cancer; CIN samples were not used) |
| GSE192897 | GEO (GSE192897) | GPL27956 (NanoString PanCancer IO 360) | 14 healthy + 12 non-recurrent primary cancers (full series: 48 = 14 healthy + 12 non-recurrent + 14 local recurrence + 8 distant recurrence; recurrent samples were not used) |

* Discovery cohort after ComBat integration: **78 samples = 38 normal + 40 cancer** (`data/merge.normalize.txt`).
* **External validation (404 samples)**: 306 TCGA-CESC tumors + 3 TCGA adjacent-normal + 95 GTEx normal cervix, all downloaded via **UCSC Xena**.
* HPV status (GSE192897 sample-level metadata): all 12 retained cancers are HPV-positive (8 × HPV16, 3 × HPV18, 1 × HPV52); healthy controls and GSE63514 have no HPV typing data.

## 4. Files shipped in this repository

* `data/merge.normalize.txt` — the integrated, ComBat-corrected 78-sample expression matrix.
* `gene_lists/interGenes.txt` — the 7-gene DEG ∩ WGCNA ∩ NETs intersection.
* `gene_lists/interFeatureGenes.txt` — the final 3-gene signature (ARG1, GZMB, EPCAM).
* `gene_lists/module_turquoise.txt` — the 148-gene WGCNA turquoise module.
* `gene_lists/NETs.txt` — NETs-related gene list derived from a GeneCards "neutrophil extracellular trap" search (the raw GeneCards export is **not** redistributed for licensing reasons).
* `reference_data/immune.gmt` — 28 immune-cell signature gene sets used for ssGSEA.
* `reference_data/geneREF.txt` — hg38 gene coordinates (chromosome, start, end) used by the circle/Manhattan plots.

## 5. Not included (must be obtained separately)

* **`data/GTExTCGA.merge.txt`** (162 MB) — external validation matrix; exceeds GitHub file limits. Build it from UCSC Xena (TCGA-CESC + GTEx cervical normal) with columns named `<source>_<sample>_<Type>` where `Type` ∈ {`Cancer`, `Normal`}.
* Raw GEO expression matrices (`geneMatrix.txt` per dataset) — download from GEO.
* MSigDB `c2.cp.kegg` gene-set file (for script 17) — download from MSigDB and save as `reference_data/c2.cp.kegg.Hs.symbols.gmt`.
* DSigDB drug-target table (for script 20) — download from the Enrichr/DSigDB resource and save the first two columns (term, gene) as `data/DSigDB_All_detailed.txt`.
* `GeneCards-SearchResults.csv` — not redistributable.
* Figure source images (TIFF/AI) — figure assembly was done in vector-graphics (AI) software and is not part of this code repository.

## 6. Honesty / provenance disclosure

The scripts in `scripts/` were **rebuilt from the original analysis working files** for publication. In doing so:

1. Obvious typographical damage from manual editing was repaired by evident intent (e.g., `pheapheatmap(` → `pheatmap(`, `def()`/`dv.off()` → `dev.off()`, `pd0000f(` → `pdf(`, a duplicated `corrplot(cor1` call, a broken `grepl(` filter in the immune-correlation script). Each repair is marked `FIXED` or `REPAIRED` in the script header or inline comment.
2. **Two methodological corrections were applied** and are documented in the response to the reviewer and in the revised manuscript:
   * **Script 04 (DEG):** the original design matrix ordered factor levels alphabetically (Cancer, Normal) but renamed the columns in the opposite order, so the "Cancer − Normal" contrast actually computed Normal − Cancer and all logFC values (and volcano-plot Up/Down labels) had reversed signs. The fix sets the factor levels explicitly. Biological conclusions are unchanged (ARG1 down-regulated in tumors; the other six genes up-regulated).
   * **Script 14 (external validation):** the original script refit the logistic model inside the validation set (AUC = 0.994 is an in-cohort refitted estimate). The corrected primary analysis applies the training-cohort coefficients unchanged (locked model: intercept 0.8263, ARG1 −1.0028, GZMB 0.3379, EPCAM 0.1328) → **AUC = 0.968 (95% CI 0.940–0.989)**; the refit is retained only as a clearly labelled sensitivity analysis.
3. **Reproducibility:** `set.seed(12345)` was added to the LASSO and SVM-RFE scripts (it was commented out in the working files); RF (12345) and consensus clustering (123456) already had seeds.
4. Per-gene ROC curves use pROC's automatic direction handling; on the external cohort ARG1's raw orientation is 0.181 and is reversed to 0.819 because ARG1 is lower in tumors.
5. Correlation analyses (script 19c) are computed **within the 40 tumor samples only**, matching the original analysis and the multiple-testing table (Supplementary Table S3B).
6. The K range reported in the manuscript (K = 2–12) is reproduced with `maxK = 12` in script 15; the surviving working file had `maxK = 9`, and the K = 3 solution (C1 = 20, C2 = 12, C3 = 8) matches the published cluster assignments.

## 7. Key quantitative results reproduced by this pipeline

* Cohort: 78 samples (38 normal, 40 cancer); external validation: 404 samples (306 + 3 + 95).
* DEGs: 123 (|log2FC| > 1, BH adj.P < 0.05).
* WGCNA: 634 genes after SD > 0.5 filter; beta = 8; turquoise module = 148 genes, correlated with tumor trait (r = 0.37, P = 8e-4).
* Signature: ARG1, GZMB, EPCAM (LASSO ∩ SVM-RFE ∩ RF).
* Training AUC = 0.864 (95% CI 0.770–0.934); Youden cutoff 0.722 → Sensitivity 65.0%, Specificity 92.1%.
* External validation, locked coefficients: AUC = 0.968 (95% CI 0.940–0.989); Youden 0.812 → 96.7% / 96.9%. Refit (sensitivity): 0.994 (0.986–1.000).
* Single-gene external AUCs: ARG1 0.819, GZMB 0.959, EPCAM 0.981.
* Subtypes (K = 3): C1 = 20, C2 = 12, C3 = 8.
* Gene × immune correlations within tumors: 154 tests, 48 with nominal P < 0.05, 22 with BH q < 0.05 (Supplementary Table S3B).

## 8. License

MIT — see [LICENSE](LICENSE).
