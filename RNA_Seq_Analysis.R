################################################################################
# Project: Transcriptomic Analysis of Alzheimer's Disease (NLGF Mouse Model)
# Model: Multi-factor Design with Age * Genotype Interaction
# Author: Daniele Brambilla
# Purpose: Identify genes showing differential expression trajectories over time
################################################################################

# --- 1. LIBRARIES & SETUP ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
libraries <- c("DESeq2", "ggplot2", "pheatmap", "GEOquery", "apeglm")
lapply(libraries, function(x) {
  if (!requireNamespace(x, quietly = TRUE)) BiocManager::install(x)
  library(x, character.only = TRUE)
})

# --- 2. DATA IMPORT & CLEANING ---
counts <- read.table("GSE290305_AD_mouse_RNAseq-rawGeneCounts_table.txt",
                     header = TRUE, row.names = 1)

gse <- getGEO(filename = "GSE290305_series_matrix.txt")
coldata_raw <- pData(gse)

# Harmonize sample IDs
new_names <- coldata_raw$title
new_names <- gsub("Hippocampus_", "", new_names)
new_names <- gsub("_RNAseq", "", new_names)
new_names <- gsub("Female", "F", new_names)
new_names <- gsub("Male", "M", new_names)
new_names <- gsub("W([0-9]+)", "T\\1", new_names)
new_names <- gsub("(.*)_(T[0-9]+)_(.*)", "\\2_\\1_\\3", new_names)
rownames(coldata_raw) <- new_names

# Align counts and metadata
common_samples <- intersect(colnames(counts), rownames(coldata_raw))
counts <- counts[, common_samples]
coldata_raw <- coldata_raw[common_samples, ]

# Final experimental design dataframe
coldata_final <- data.frame(
  sample_id = colnames(counts),
  genotype  = factor(ifelse(grepl("WT", colnames(counts)), "WT", "NLGF")),
  age       = factor(gsub("(T[0-9]+)_.*", "\\1", colnames(counts)),
                     levels = c("T3", "T8", "T24")),
  sex       = factor(ifelse(grepl("_F", colnames(counts)), "Female", "Male")),
  row.names = colnames(counts)
)

coldata_final$genotype <- relevel(coldata_final$genotype, ref = "WT")
stopifnot(all(colnames(counts) == rownames(coldata_final)))

# --- 3. DIFFERENTIAL EXPRESSION: INTERACTION MODEL ---

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata_final,
  design    = ~ age + sex + genotype + age:genotype
)

# --- 3A. LIKELIHOOD RATIO TEST (LRT) ---
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ age + sex + genotype)
res_lrt <- results(dds_lrt)
res_lrt <- res_lrt[order(res_lrt$padj), ]

# NEW: Filter out NAs before selecting significant genes
# This prevents the "logical subscript contains NAs" error
sig_lrt_results <- res_lrt[!is.na(res_lrt$padj) & res_lrt$padj < 0.05, ]
sig_lrt_genes <- rownames(sig_lrt_results)

# --- 3B. WALD TEST (WHERE / WHEN EFFECT OCCURS)
# Used to interpret genotype effects at specific ages (e.g. T24)

dds_wald <- DESeq(dds, test = "Wald")
resultsNames(dds_wald)

# Effect of NLGF vs WT specifically at T24
res_T24 <- results(
  dds_wald,
  contrast = list(
    c("genotype_NLGF_vs_WT", "ageT24.genotypeNLGF")
  )
)

res_T24 <- res_T24[order(res_T24$padj), ]

# --- 4. VARIANCE STABILIZATION & PCA ---

vsd <- vst(dds_wald, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = c("genotype", "age"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = genotype, shape = age)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCA: Age × Genotype Interaction",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  )

# --- 5. HEATMAP OF INTERACTION GENES (LRT) ---
# Select the top 30 most significant interaction genes
top_lrt_genes <- head(sig_lrt_genes, 30)

# Check if we have enough genes, then plot
if(length(top_lrt_genes) > 0) {
  pheatmap(
    assay(vsd)[top_lrt_genes, ],
    scale = "row",
    clustering_distance_rows = "correlation",
    annotation_col = coldata_final[, c("genotype", "age", "sex")],
    main = "Genes with Significant Age × Genotype Interaction",
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
  )
} else {
  print("No significant interaction genes found with padj < 0.05")
}
# --- 6. EXPORT RESULTS ---

write.csv(as.data.frame(res_lrt),
          "AD_LRT_Age_Genotype_Interaction_AllGenes.csv")

write.csv(as.data.frame(res_T24),
          "AD_Wald_NLGF_vs_WT_T24.csv")
