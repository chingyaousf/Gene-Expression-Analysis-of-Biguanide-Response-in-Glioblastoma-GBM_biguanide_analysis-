
# -----------------------------------------------
# Annotate DESeq2 Results with Gene Symbols (Full Table)
# -----------------------------------------------

# Load required libraries
library(org.Hs.eg.db)
library(AnnotationDbi)

# Load DESeq2 results (full)
res_df <- read.csv("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/deseq2_results.csv")

# Extract Ensembl IDs without version
res_df$ENSEMBL_clean <- sub("\\..*", "", res_df$Gene)

# Map Ensembl IDs to gene symbols
symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(res_df$ENSEMBL_clean),
  keytype = "ENSEMBL",
  columns = c("SYMBOL")
)

# Remove duplicated mappings (keep one symbol per Ensembl)
symbol_map <- symbol_map[!duplicated(symbol_map$ENSEMBL), ]

# Merge symbol into DESeq2 results
res_df$SYMBOL <- symbol_map$SYMBOL[match(res_df$ENSEMBL_clean, symbol_map$ENSEMBL)]

# Optional: move SYMBOL column next to Gene column
res_df <- res_df[, c("Gene", "SYMBOL", setdiff(names(res_df), c("Gene", "SYMBOL")))]

# Save updated results
write.csv(res_df, "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/deseq2_results_with_symbols.csv", row.names = FALSE)





# -----------------------------------------------
# Heatmap of Top 50 Differentially Expressed Genes (Using Gene Symbols)
# -----------------------------------------------

# Load required libraries
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(tibble)
library(org.Hs.eg.db)
library(AnnotationDbi)

# -----------------------------------------------
# Load normalized counts and DESeq2 results
# -----------------------------------------------
norm_counts <- read.csv("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/normalized_counts.csv", row.names = 1)
#res_df <- read.csv("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/deseq2_results.csv")
res_df <- readRDS("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/res_df_annotated.rds")

# -----------------------------------------------
# Select top 50 DEGs based on adjusted p-value
# -----------------------------------------------
top_degs <- res_df[!is.na(res_df$padj), ]
top_degs <- top_degs[order(top_degs$padj), ]
top50_genes <- head(top_degs$Gene, 50)  # Ensembl IDs with version

# -----------------------------------------------
# Clean Ensembl IDs (remove version numbers)
# -----------------------------------------------
top50_clean <- sub("\\..*", "", top50_genes)

# Map Ensembl to gene symbols
symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = top50_clean,
  keytype = "ENSEMBL",
  columns = c("SYMBOL")
)

# Remove duplicates and NA
symbol_map <- symbol_map[!duplicated(symbol_map$ENSEMBL) & !is.na(symbol_map$SYMBOL), ]

# Subset normalized counts using original Ensembl IDs
heatmap_mat <- norm_counts[top50_genes, ]

# Scale by row (z-score)
heatmap_scaled <- t(scale(t(heatmap_mat)))

# -----------------------------------------------
# Replace Ensembl IDs with gene symbols in heatmap rows
# -----------------------------------------------
rownames(heatmap_scaled) <- symbol_map$SYMBOL[match(sub("\\..*", "", rownames(heatmap_scaled)), symbol_map$ENSEMBL)]

# -----------------------------------------------
# Create sample annotation (Sensitive vs Resistant)
# -----------------------------------------------
sample_info <- data.frame(
  Condition = c("Sensitive", "Sensitive", "Sensitive", "Resistant", "Resistant", "Sensitive")
)
rownames(sample_info) <- colnames(heatmap_scaled)

ann_colors <- list(
  Condition = c(Sensitive = "#1B9E77", Resistant = "#D95F02")
)

# -----------------------------------------------
# Plot and Save Heatmap of Top 50 DEGs (with Gene Symbols)
# -----------------------------------------------
png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/heatmap_top50_DEGs_02.png", 
    width = 1000, height = 1200, res = 150)

pheatmap(
  heatmap_scaled,
  annotation_col = sample_info,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 12,
  main = "Top 50 DEGs (Resistant vs Sensitive)",
  color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)
)

dev.off()




# -----------------------------------------------
# Heatmap of Epigenetic-Related Genes from GSEA (SYMBOL-based)
# -----------------------------------------------

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Step 1: Load normalized counts and GSEA results
norm_counts <- read.csv("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/normalized_counts.csv", row.names = 1)
gsea_res <- read.csv("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/gsea_epigenetic_results.csv")

# Step 2: Convert rownames (ENSEMBL) → SYMBOL
ensembl_ids <- sub("\\..*", "", rownames(norm_counts))  # Remove version numbers
symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ensembl_ids,
  keytype = "ENSEMBL",
  columns = "SYMBOL"
)

# Merge SYMBOLs into norm_counts
norm_counts$ENSEMBL <- ensembl_ids
norm_counts <- left_join(norm_counts, symbol_map, by = c("ENSEMBL" = "ENSEMBL"))
norm_counts <- norm_counts[!duplicated(norm_counts$SYMBOL) & !is.na(norm_counts$SYMBOL), ]
rownames(norm_counts) <- norm_counts$SYMBOL
norm_counts <- norm_counts[, c("PDX22", "PDX43", "PDX59", "G83", "PDX10", "HK281")]

# Step 3: Get gene symbols from top GSEA epigenetic pathways
top_paths <- gsea_res %>% arrange(p.adjust) %>% head(10)
epigenetic_genes <- unique(unlist(strsplit(top_paths$core_enrichment, "/")))

# Step 4: Filter normalized counts
epi_heatmap_mat <- norm_counts[rownames(norm_counts) %in% epigenetic_genes, ]

# Step 5: Z-score transform
epi_heatmap_scaled <- t(scale(t(epi_heatmap_mat)))

# Step 6: Add sample annotation
sample_info <- data.frame(
  Condition = c("Sensitive", "Sensitive", "Sensitive", "Resistant", "Resistant", "Sensitive")
)
rownames(sample_info) <- colnames(epi_heatmap_scaled)
ann_colors <- list(Condition = c(Sensitive = "#1B9E77", Resistant = "#D95F02"))

# Step 7: Plot heatmap
png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/heatmap_epigenetic_genes.png",
    width = 1000, height = 1200, res = 150)

pheatmap(
  epi_heatmap_scaled,
  annotation_col = sample_info,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  fontsize_row = 9,
  fontsize_col = 12,
  main = "Epigenetic Genes (Top GSEA Pathways)",
  color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)
)

dev.off()



# -----------------------------------------------
# Heatmap of Differentially Expressed Transcription Factors
# -----------------------------------------------

library(dplyr)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tibble)
library(dorothea)

# Step 1: Load normalized counts and DE results
norm_counts <- read.csv("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/normalized_counts.csv", row.names = 1)
res_df <- readRDS("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/res_df_annotated.rds")

# Step 2: Map normalized counts from ENSEMBL → SYMBOL
ensembl_ids <- sub("\\..*", "", rownames(norm_counts))
symbol_map <- AnnotationDbi::select(org.Hs.eg.db, keys = ensembl_ids, keytype = "ENSEMBL", columns = "SYMBOL")

# Attach SYMBOL and clean
norm_counts$ENSEMBL <- ensembl_ids
norm_counts <- left_join(norm_counts, symbol_map, by = c("ENSEMBL" = "ENSEMBL"))
norm_counts <- norm_counts[!duplicated(norm_counts$SYMBOL) & !is.na(norm_counts$SYMBOL), ]
rownames(norm_counts) <- norm_counts$SYMBOL
norm_counts <- norm_counts[, c("PDX22", "PDX43", "PDX59", "G83", "PDX10", "HK281")]

# Step 3: Get curated TF list from dorothea
data(dorothea_hs)
tf_list <- unique(dorothea_hs$tf)

# Step 4: Get significant DEGs that are TFs
deg_tfs <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1 & SYMBOL %in% tf_list)

# Step 5: Match those TFs to normalized counts
tf_genes <- deg_tfs$SYMBOL
tf_heatmap_mat <- norm_counts[rownames(norm_counts) %in% tf_genes, ]

# Z-score scale
tf_heatmap_scaled <- t(scale(t(tf_heatmap_mat)))

# Step 6: Add sample annotation
sample_info <- data.frame(
  Condition = c("Sensitive", "Sensitive", "Sensitive", "Resistant", "Resistant", "Sensitive")
)
rownames(sample_info) <- colnames(tf_heatmap_scaled)
ann_colors <- list(Condition = c(Sensitive = "#1B9E77", Resistant = "#D95F02"))

# Step 7: Plot the heatmap
png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/heatmap_tf_DEGs.png", 
    width = 1000, height = 1200, res = 150)

pheatmap(
  tf_heatmap_scaled,
  annotation_col = sample_info,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  fontsize_row = 9,
  fontsize_col = 12,
  main = "Differentially Expressed Transcription Factors",
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = "RdBu")))(100)
)

dev.off()



##  Heatmap of Epigenetic-Related Transcription Factors
# ------------------------------------------------------
# Step 1: Load normalized counts and epigenetic GSEA results
# ------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(dorothea)
library(pheatmap)
library(RColorBrewer)
library(tibble)

# Load data
norm_counts <- read.csv("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/normalized_counts.csv", row.names = 1)
gsea_epigenetic <- read.csv("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/gsea_epigenetic_results.csv")

# ------------------------------------------------------
# Step 2: Extract unique genes from `core_enrichment`
# ------------------------------------------------------
epigenetic_genes <- unique(unlist(strsplit(gsea_epigenetic$core_enrichment, "/")))

# ------------------------------------------------------
# Step 3: Load DoRothEA transcription factor list
# ------------------------------------------------------
data(dorothea_hs, package = "dorothea")
tfs_all <- unique(dorothea_hs$target)

# Intersect TFs with epigenetic pathway genes
epigenetic_tfs <- intersect(epigenetic_genes, tfs_all)

# ------------------------------------------------------
# Step 4: Map SYMBOL → ENSEMBL
# ------------------------------------------------------
symbol_to_ens <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = epigenetic_tfs,
  keytype = "SYMBOL",
  columns = "ENSEMBL"
)
symbol_to_ens <- symbol_to_ens[!duplicated(symbol_to_ens$SYMBOL), ]

# ------------------------------------------------------
# Step 5: Clean ENSEMBL rownames in normalized matrix
# ------------------------------------------------------
ensembl_clean <- sub("\\..*", "", rownames(norm_counts))
norm_counts$ensembl_clean <- ensembl_clean

# Collapse duplicates by mean (can change to median or first)
norm_counts_dedup <- aggregate(. ~ ensembl_clean, data = norm_counts, FUN = mean)

# Set cleaned rownames
rownames(norm_counts_dedup) <- norm_counts_dedup$ensembl_clean
norm_counts_dedup$ensembl_clean <- NULL

# ------------------------------------------------------
# Step 6: Subset and scale TF expression matrix
# ------------------------------------------------------
matched_ensembl <- intersect(symbol_to_ens$ENSEMBL, rownames(norm_counts_dedup))
epi_tf_mat <- norm_counts_dedup[matched_ensembl, ]
epi_tf_scaled <- t(scale(t(epi_tf_mat)))

# Rename rows to gene symbols
rownames(epi_tf_scaled) <- symbol_to_ens$SYMBOL[match(rownames(epi_tf_scaled), symbol_to_ens$ENSEMBL)]

# ------------------------------------------------------
# Step 7: Annotation metadata for columns
# ------------------------------------------------------
sample_info <- data.frame(
  Condition = c("Sensitive", "Sensitive", "Sensitive", "Resistant", "Resistant", "Sensitive")
)
rownames(sample_info) <- colnames(epi_tf_scaled)

ann_colors <- list(
  Condition = c(Sensitive = "#1B9E77", Resistant = "#D95F02")
)

# ------------------------------------------------------
# Step 8: Plot and Save Heatmap
# ------------------------------------------------------
png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/heatmap_epigenetic_TFs.png",
    width = 1000, height = 1200, res = 150)

pheatmap(
  epi_tf_scaled,
  annotation_col = sample_info,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 12,
  main = "Epigenetic-Related Transcription Factors (All)",
  color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)
)

dev.off()





# ------------------------------------------------------
# Full Script: Heatmaps and Volcano Plots for GBM Gene Sets
# ------------------------------------------------------

# Load required libraries
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(tibble)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)
library(dplyr)
library(dorothea)
library(clusterProfiler)
library(ggplot2)
library(ggrepel)

# ------------------------------------------------------
# Volcano Plot Function (Reusable)
# ------------------------------------------------------
plot_volcano <- function(res_df, highlight_genes = NULL, title = "Volcano Plot", filename = "volcano.png") {
  res_df <- res_df %>%
    mutate(
      log10padj = -log10(padj),
      Significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"),
      Label = ifelse(SYMBOL %in% highlight_genes, SYMBOL, NA)
    )

  p <- ggplot(res_df, aes(x = log2FoldChange, y = log10padj)) +
    geom_point(aes(color = Significant), alpha = 0.4) +
    geom_point(data = subset(res_df, SYMBOL %in% highlight_genes), color = "red", size = 0.5) +
    geom_text_repel(aes(label = Label), max.overlaps = 30) +
    scale_color_manual(values = c("No" = "gray70", "Yes" = "#1B9E77")) +
    labs(title = title, x = "log2 Fold Change", y = "-log10 Adjusted P-value") +
    theme_minimal()

  ggsave(filename, plot = p, width = 7, height = 6, dpi = 300)
}

# ------------------------------------------------------
# Load data
# ------------------------------------------------------
norm_counts <- read.csv("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/normalized_counts.csv", row.names = 1)
res_df <- readRDS("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/res_df_annotated.rds")
gsea_epigenetic <- read.csv("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/gsea_epigenetic_results.csv")

# ------------------------------------------------------
# 1. Top 50 DEGs Volcano Plot
# ------------------------------------------------------
top_degs <- res_df[!is.na(res_df$padj), ] %>% arrange(padj) %>% head(50)
top50_genes <- sub("\\..*", "", top_degs$Gene)
symbol_map_top50 <- AnnotationDbi::select(org.Hs.eg.db, keys = top50_genes, keytype = "ENSEMBL", columns = "SYMBOL")
symbol_map_top50 <- symbol_map_top50[!duplicated(symbol_map_top50$ENSEMBL) & !is.na(symbol_map_top50$SYMBOL), ]

plot_volcano(
  res_df = res_df,
  highlight_genes = symbol_map_top50$SYMBOL,
  title = "Top 50 DEGs",
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/volcano_top50_DEGs.png"
)

# ------------------------------------------------------
# 2. Epigenetic Genes Volcano Plot
# ------------------------------------------------------
top_paths <- gsea_epigenetic %>% arrange(p.adjust) %>% head(10)
epigenetic_genes <- unique(unlist(strsplit(top_paths$core_enrichment, "/")))

plot_volcano(
  res_df = res_df,
  highlight_genes = epigenetic_genes,
  title = "Epigenetic Genes from GSEA",
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/volcano_epigenetic_genes.png"
)

# ------------------------------------------------------
# 3. Differentially Expressed Transcription Factors Volcano Plot
# ------------------------------------------------------
data(dorothea_hs)
tf_list <- unique(dorothea_hs$tf)
deg_tfs <- res_df %>% filter(padj < 0.05 & abs(log2FoldChange) > 1 & SYMBOL %in% tf_list)

plot_volcano(
  res_df = res_df,
  highlight_genes = deg_tfs$SYMBOL,
  title = "Differentially Expressed TFs",
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/volcano_tf_DEGs.png"
)

# ------------------------------------------------------
# 4. Epigenetic-Related Transcription Factors Volcano Plot
# ------------------------------------------------------
tfs_all <- unique(dorothea_hs$target)
epigenetic_tfs <- intersect(epigenetic_genes, tfs_all)

plot_volcano(
  res_df = res_df,
  highlight_genes = epigenetic_tfs,
  title = "Epigenetic-Related TFs",
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/volcano_epigenetic_TFs.png"
)



# using EnhancedVolcano 
# -----------------------------------------------
# Load Required Libraries
# -----------------------------------------------
library(EnhancedVolcano)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dorothea)

# -----------------------------------------------
# Load Data
# -----------------------------------------------
norm_counts <- read.csv("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/normalized_counts.csv", row.names = 1)
res_df <- readRDS("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/res_df_annotated.rds")
gsea_epigenetic <- read.csv("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/gsea_epigenetic_results.csv")

# -----------------------------------------------
# 1. Volcano Plot: Top 50 DEGs
# -----------------------------------------------
top_degs <- res_df[!is.na(res_df$padj), ] %>% arrange(padj) %>% head(50)
top50_genes <- sub("\\..*", "", top_degs$Gene)

symbol_map_top50 <- AnnotationDbi::select(org.Hs.eg.db, keys = top50_genes, keytype = "ENSEMBL", columns = "SYMBOL")
symbol_map_top50 <- symbol_map_top50[!duplicated(symbol_map_top50$ENSEMBL) & !is.na(symbol_map_top50$SYMBOL), ]

EnhancedVolcano(res_df,
  lab = res_df$SYMBOL,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = symbol_map_top50$SYMBOL,
  title = 'Top 50 DEGs',
  subtitle = 'Resistant vs Sensitive',
  pCutoff = 0.05,
  FCcutoff = 1,
  colAlpha = 0.5,
  pointSize = 2,
  labSize = 4,
  max.overlaps = Inf,
  drawConnectors = TRUE,
  widthConnectors = 0.3,
  arrowheads = FALSE 
)
ggsave("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/volcano_top50_DEGs_enhanced_02.png", width = 15, height = 10, dpi = 300)

# -----------------------------------------------
# 2. Volcano Plot: Epigenetic Genes from GSEA
# -----------------------------------------------
top_paths <- gsea_epigenetic %>% arrange(p.adjust) %>% head(10)
epigenetic_genes <- unique(unlist(strsplit(top_paths$core_enrichment, "/")))

EnhancedVolcano(res_df,
  lab = res_df$SYMBOL,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = epigenetic_genes,
  title = 'Epigenetic Genes from GSEA',
  subtitle = 'Top Pathways (Reactome)',
  pCutoff = 0.05,
  FCcutoff = 1,
  colAlpha = 0.5,
  pointSize = 2,
  labSize = 4,
  max.overlaps = Inf,
  drawConnectors = TRUE,
  widthConnectors = 0.3,
  arrowheads = FALSE 
)
ggsave("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/volcano_epigenetic_genes_enhanced_02.png", width = 15, height = 10, dpi = 300)




# -----------------------------------------------
# 2. Volcano Plot: Epigenetic Genes from GSEA (Significant Only)
# -----------------------------------------------

library(EnhancedVolcano)
library(dplyr)

# Load data
res_df <- readRDS("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/res_df_annotated.rds")
gsea_epigenetic <- read.csv("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/gsea_epigenetic_results.csv")

# Extract epigenetic genes from top 10 GSEA pathways
top_paths <- gsea_epigenetic %>% arrange(p.adjust) %>% head(10)
epigenetic_genes <- unique(unlist(strsplit(top_paths$core_enrichment, "/")))

# Filter for epigenetic genes that are significantly differentially expressed
epigenetic_significant <- res_df %>%
  filter(SYMBOL %in% epigenetic_genes & padj < 0.05 & abs(log2FoldChange) > 1) %>%
  pull(SYMBOL)

# Plot using EnhancedVolcano
EnhancedVolcano(res_df,
  lab = res_df$SYMBOL,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = epigenetic_significant,
  title = 'Epigenetic Genes from GSEA',
  subtitle = 'Top Pathways (Reactome)',
  pCutoff = 0.05,
  FCcutoff = 1,
  colAlpha = 0.5,
  pointSize = 2,
  labSize = 4,
  max.overlaps = Inf,
  drawConnectors = TRUE,
  widthConnectors = 0.3,
  arrowheads = FALSE
)

# Save plot
ggsave("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/volcano_epigenetic_genes_enhanced_significant_02.png",
       width = 15, height = 10, dpi = 300)


# -----------------------------------------------
# 3. Volcano Plot: Differentially Expressed TFs
# -----------------------------------------------
data(dorothea_hs)
tf_list <- unique(dorothea_hs$tf)

deg_tfs <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1 & SYMBOL %in% tf_list)

EnhancedVolcano(res_df,
  lab = res_df$SYMBOL,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = deg_tfs$SYMBOL,
  title = 'Differentially Expressed TFs',
  pCutoff = 0.05,
  FCcutoff = 1,
  colAlpha = 0.5,
  pointSize = 2,
  labSize = 4,
  max.overlaps = Inf,
  drawConnectors = TRUE,
  widthConnectors = 0.3,
  arrowheads = FALSE
)
ggsave("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/volcano_tf_DEGs_enhanced_02.png", width = 15, height = 10, dpi = 300)

# -----------------------------------------------
# 4. Volcano Plot: Epigenetic-Related TFs
# -----------------------------------------------
tfs_all <- unique(dorothea_hs$target)
epigenetic_tfs <- intersect(epigenetic_genes, tfs_all)

EnhancedVolcano(res_df,
  lab = res_df$SYMBOL,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = epigenetic_tfs,
  title = 'Epigenetic-Related TFs',
  pCutoff = 0.05,
  FCcutoff = 1,
  colAlpha = 0.5,
    pointSize = 2,
  labSize = 4,
  max.overlaps = Inf,
  drawConnectors = TRUE,
  widthConnectors = 0.3,
  arrowheads = FALSE
)
ggsave("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/volcano_epigenetic_TFs_enhanced_02.png", width = 15, height = 10, dpi = 300)
