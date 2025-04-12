
# ================================================
# Step 1: Load and Prepare Count Matrix for Analysis
# ================================================

# Load raw count matrix from featureCounts output
counts_raw <- read.delim("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/counts_paired.txt", 
                         comment.char = "#",        # Ignore lines starting with '#' (featureCounts metadata)
                         header = TRUE,             # First row is header
                         row.names = 1,             # First column is gene ID
                         check.names = FALSE)       # Keep sample names as-is (don't auto-fix them)

# Remove annotation columns: Chr, Start, End, Strand, Length
# Keep only columns that contain actual count data
counts <- counts_raw[, 6:ncol(counts_raw)]

# Rename columns to short, readable sample names
colnames(counts) <- c("PDX22", "PDX43", "PDX59", "G83", "PDX10", "HK281")

# Preview the cleaned count matrix
print(head(counts))

# Create sample metadata indicating experimental groups
sample_info <- data.frame(
  sample = colnames(counts),
  condition = c("Sensitive", "Sensitive", "Sensitive", "Resistant", "Resistant", "Sensitive")  # Group assignments
)

# Set row names to match sample names
rownames(sample_info) <- sample_info$sample

# Preview metadata table
print(sample_info)

# Save cleaned count matrix and metadata for downstream analysis
write.csv(counts, "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/cleaned_counts.csv")
write.csv(sample_info, "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/sample_metadata.csv", row.names = FALSE)




# -----------------------------------------------
# Step 2: Differential Expression Analysis (DESeq2)
# -----------------------------------------------

# Load required libraries
library(DESeq2)
library(tibble)

# Read in cleaned count matrix and metadata
counts <- read.csv("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/cleaned_counts.csv", row.names = 1)
metadata <- read.csv("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/sample_metadata.csv", row.names = 1)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = metadata,
                              design = ~ condition)

# Run DESeq2
dds <- DESeq(dds)

# -----------------------------------------------
# Extract and Save Normalized Counts
# -----------------------------------------------
normalized_counts <- counts(dds, normalized = TRUE)
normalized_counts_rounded <- round(normalized_counts, 2)
write.csv(normalized_counts_rounded, "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/normalized_counts.csv")

# -----------------------------------------------
# Get Differential Expression Results
# -----------------------------------------------
# Get results for comparison: Resistant vs Sensitive
res <- results(dds, contrast = c("condition", "Resistant", "Sensitive"))

# Save the DESeq2 results object (res) for later use (e.g., for plotMA)
saveRDS(res, file = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/res_object.rds")

# Order by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Convert to data frame and add gene names
res_df <- as.data.frame(res_ordered)
res_df <- rownames_to_column(res_df, var = "Gene")

# Save results to CSV
write.csv(res_df, "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/deseq2_results.csv", row.names = FALSE)

# Optional: Summary output
summary(res)

# Save the DESeq2 dataset object for future use (e.g., QC plots, VST)
saveRDS(dds, file = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/dds_object.rds")


# -----------------------------------------------
# Step 3: QC Plots After DESeq2 Normalization (with PNG Output)
# -----------------------------------------------

# Load required libraries
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

# Create output directory
output_dir <- "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/qc_plots"
dir.create(output_dir, showWarnings = FALSE)

# Load DESeq2 object saved from Step 2
dds <- readRDS("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/dds_object.rds")

# VST Transformation
vsd <- vst(dds, blind = FALSE)

# Get results for comparison: Resistant vs Sensitive
res <- results(dds, contrast = c("condition", "Resistant", "Sensitive"))

# 1. Sample Distance Heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)

png(file.path(output_dir, "sample_distance_heatmap.png"), width = 1000, height = 1000)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colorRampPalette(rev(brewer.pal(9, "Greens")))(255),
         main = "Sample Distance Heatmap",
         fontsize_row = 18,    # Y-axis label size
         fontsize_col = 18,    # X-axis label size
         fontsize = 19)        # General font size (e.g. legend)
dev.off()


# 2. PCA Plot (with sample labels and bigger fonts)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$SampleName <- rownames(pcaData)

p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, label = SampleName)) +
  geom_point(size = 5) +
  geom_text(vjust = -1.2, size = 4) +  # Add sample names
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: Resistant vs Sensitive") +
  theme_bw(base_size = 18) +  # Bigger fonts
  theme(legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 20))

# Save high-res PCA plot
png(file.path(output_dir, "PCA_plot.png"), width = 1400, height = 1200, res = 150)
print(p)
dev.off()


# 3. MA Plot
png(file.path(output_dir, "MA_plot.png"), width = 800, height = 800)
plotMA(res, main = "MA Plot: Resistant vs Sensitive", ylim = c(-5, 5))
dev.off()

# 3. MA Plot with larger labels
png(file.path(output_dir, "MA_plot_big_labels.png"), width = 1000, height = 1000, res = 150)

# Adjust base R plot parameters:
par(
  cex = 1.2,        # Overall scaling of text and symbols
  cex.axis = 1.2,   # Axis tick label size
  cex.lab = 1.2,    # Axis label size
  cex.main = 1.4    # Main title size
)

plotMA(
  res,
  main = "MA Plot: Resistant vs Sensitive",
  ylim = c(-5, 5),
  cex = 0.3          # Point size
)
dev.off()

# 4. Dispersion Estimates
png(file.path(output_dir, "dispersion_plot.png"), width = 800, height = 800)
plotDispEsts(dds, main = "Dispersion Estimates")
dev.off()

# 4. Dispersion Estimates (with larger labels)
png(file.path(output_dir, "dispersion_plot.png"), width = 1000, height = 1000, res = 150)

# Set font sizes for base R plot
par(
  cex = 1.2,        # Overall scaling
  cex.axis = 1.2,   # Axis tick labels
  cex.lab = 1.4,    # Axis titles
  cex.main = 1.6    # Plot title
)

plotDispEsts(dds, main = "Dispersion Estimates")

dev.off()




# -----------------------------------------------
# Step 4: Define Upregulated and Downregulated DEGs
# -----------------------------------------------

# Load DESeq2 result object
res <- readRDS("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/res_object.rds")

# Convert to data frame and add gene names
res_df <- as.data.frame(res)
res_df <- tibble::rownames_to_column(res_df, var = "Gene")

# Now define up/downregulated DEGs
up_genes <- subset(res_df, log2FoldChange > 1 & padj < 0.05)$Gene
down_genes <- subset(res_df, log2FoldChange < -1 & padj < 0.05)$Gene

# Strip version numbers from Ensembl IDs
up_genes_clean <- sub("\\..*", "", up_genes)
down_genes_clean <- sub("\\..*", "", down_genes)

head(up_genes_clean)
head(down_genes_clean)

# -----------------------------------------------
# Step 5: GO Enrichment Analysis (clusterProfiler)
# -----------------------------------------------

# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)

# -----------------------------------------------
# Step 5A: Convert cleaned Ensembl IDs to Entrez IDs
# -----------------------------------------------
gene_up_entrez <- bitr(up_genes_clean, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
gene_down_entrez <- bitr(down_genes_clean, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

#check how many genes successfully converted
length(up_genes_clean); length(gene_up_entrez)
length(down_genes_clean); length(gene_down_entrez)


# -----------------------------------------------
# Step 5B: GO Biological Process Enrichment
# -----------------------------------------------
ego_up <- enrichGO(
  gene = gene_up_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  #qvalueCutoff = 0.1,
  readable = TRUE
)

# See how many terms are enriched
nrow(as.data.frame(ego_up))

ego_down <- enrichGO(
  gene = gene_down_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  #qvalueCutoff = 0.1,
  readable = TRUE
)

# See how many terms are enriched
nrow(as.data.frame(ego_down))


# View top 10 enriched BP terms for up-regulated genes
head(ego_up, 10)

# View top 10 enriched BP terms for down-regulated genes
head(ego_down, 10)

# -----------------------------------------------
# Step 5C: Save GO Enrichment Results to CSV
# -----------------------------------------------
dir.create("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment", showWarnings = FALSE)

write.csv(as.data.frame(ego_up),
          "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/go_up.csv",
          row.names = FALSE)

write.csv(as.data.frame(ego_down),
          "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/go_down.csv",
          row.names = FALSE)

# -----------------------------------------------
# Step 5D: GO Enrichment Dotplots
# -----------------------------------------------
# Create directory if it doesn't exist
dir.create("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots", showWarnings = FALSE)

# ---- Dotplot: Upregulated Genes ----
png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/dotplot_go_up.png", width = 600, height = 400)
dotplot(ego_up, showCategory = 20, title = "GO Enrichment: Upregulated Genes")
dev.off()

# ---- Dotplot: Downregulated Genes ----
png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/dotplot_go_down.png", width = 600, height = 400)
dotplot(ego_down, showCategory = 20, title = "GO Enrichment: Downregulated Genes")
dev.off()

# ---- Barplot: Upregulated Genes ----
png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/barplot_go_up.png", width = 600, height = 400)
barplot(ego_up, showCategory = 20, title = "GO Enrichment: Upregulated Genes")
dev.off()

# ---- Barplot: Downregulated Genes ----
png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/barplot_go_down.png", width = 600, height = 400)
barplot(ego_down, showCategory = 20, title = "GO Enrichment: Downregulated Genes")
dev.off()


# -----------------------------------------------
# Step 6A: KEGG Enrichment Analysis
# -----------------------------------------------
library(clusterProfiler)

ekegg_up <- enrichKEGG(
  gene = gene_up_entrez,
  organism = "hsa",
  pvalueCutoff = 0.1
)

ekegg_down <- enrichKEGG(
  gene = gene_down_entrez,
  organism = "hsa",
  pvalueCutoff = 0.1
)

# Save KEGG results
write.csv(as.data.frame(ekegg_up), "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/kegg_up.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg_down), "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/kegg_down.csv", row.names = FALSE)

# Dotplots
png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/dotplot_kegg_up.png", width = 600, height = 400)
dotplot(ekegg_up, showCategory = 20, title = "KEGG Pathways: Upregulated Genes")
dev.off()

png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/dotplot_kegg_down.png", width = 600, height = 400)
dotplot(ekegg_down, showCategory = 20, title = "KEGG Pathways: Downregulated Genes")
dev.off()

# Barplots
png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/barplot_kegg_up.png", width = 600, height = 400)
barplot(ekegg_up, showCategory = 20, title = "KEGG Pathways: Upregulated Genes")
dev.off()

png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/barplot_kegg_down.png", width = 600, height = 400)
barplot(ekegg_down, showCategory = 20, title = "KEGG Pathways: Downregulated Genes")
dev.off()


# -----------------------------------------------
# Step 6B: Reactome Enrichment Analysis
# -----------------------------------------------
library(ReactomePA)

ereactome_up <- enrichPathway(
  gene = gene_up_entrez,
  organism = "human",
  pvalueCutoff = 0.1,
  readable = TRUE
)

ereactome_down <- enrichPathway(
  gene = gene_down_entrez,
  organism = "human",
  pvalueCutoff = 0.1,
  readable = TRUE
)

# Save Reactome results
write.csv(as.data.frame(ereactome_up), "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/reactome_up.csv", row.names = FALSE)
write.csv(as.data.frame(ereactome_down), "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/reactome_down.csv", row.names = FALSE)

# Dotplots
png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/dotplot_reactome_up.png", width = 600, height = 400)
dotplot(ereactome_up, showCategory = 20, title = "Reactome Pathways: Upregulated Genes")
dev.off()

png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/dotplot_reactome_down.png", width = 600, height = 400)
dotplot(ereactome_down, showCategory = 20, title = "Reactome Pathways: Downregulated Genes")
dev.off()

# Barplots
png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/barplot_reactome_up.png", width = 600, height = 400)
barplot(ereactome_up, showCategory = 20, title = "Reactome Pathways: Upregulated Genes")
dev.off()

png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/barplot_reactome_down.png", width = 600, height = 400)
barplot(ereactome_down, showCategory = 20, title = "Reactome Pathways: Downregulated Genes")
dev.off()


nrow(as.data.frame(ekegg_up))       # should be >0 to plot
nrow(as.data.frame(ekegg_down))
nrow(as.data.frame(ereactome_up))
nrow(as.data.frame(ereactome_down)) # already worked



# -----------------------------------------------
# Step 7: Transcription Factor (TF) Enrichment Analysis using DoRothEA + decoupleR
# -----------------------------------------------

# -----------------------------------------------
# Step 7A: Load Required Packages
# -----------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("decoupleR", quietly = TRUE)) BiocManager::install("decoupleR")
if (!requireNamespace("dorothea", quietly = TRUE)) BiocManager::install("dorothea")
BiocManager::install("viper")

library(decoupleR)
library(dorothea)
library(viper)

# -----------------------------------------------
# Step 7B: Prepare log2FC Vector (Using Gene Symbols)
# -----------------------------------------------
# Strip version from Ensembl IDs
res_df$Ensembl <- sub("\\..*", "", res_df$Gene)

# Map Ensembl IDs to gene symbols using org.Hs.eg.db
symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = res_df$Ensembl,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
)

# Merge gene symbols into res_df
res_df_annotated <- merge(res_df, symbol_map, by.x = "Ensembl", by.y = "ENSEMBL")

# Keep rows with valid SYMBOLs
res_df_annotated <- res_df_annotated[!is.na(res_df_annotated$SYMBOL), ]

# Save for reuse later (e.g., for GSEA)
saveRDS(res_df_annotated, file = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/res_df_annotated.rds")

# Load the annotated DESeq2 results
res_df_annotated <- readRDS("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/res_df_annotated.rds")

# Save as CSV file
write.csv(res_df_annotated, "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/res_df_annotated.csv", row.names = FALSE)

# Build named vector of log2FC using SYMBOLs
logFC_vector <- res_df_annotated$log2FoldChange
names(logFC_vector) <- res_df_annotated$SYMBOL

# Remove NA/Inf values
logFC_vector <- logFC_vector[is.finite(logFC_vector)]


# -----------------------------------------------
# Step 7C (Fixed): Format DoRothEA Regulons for decoupleR
# -----------------------------------------------
# Load DoRothEA and keep only A/B confidence
data(dorothea_hs, package = "dorothea")
regulons <- subset(dorothea_hs, confidence %in% c("A", "B"))

# Rename columns to match decoupleR expected format
# source = TF, target = target gene, mor = mode of regulation, likelihood = confidence
colnames(regulons)[colnames(regulons) == "tf"] <- "source"
colnames(regulons)[colnames(regulons) == "target"] <- "target"
colnames(regulons)[colnames(regulons) == "mor"] <- "mor"
colnames(regulons)[colnames(regulons) == "confidence"] <- "likelihood"


# -----------------------------------------------
# Step 7D: Run TF Enrichment Analysis with VIPER
# -----------------------------------------------
tf_activity <- decoupleR::run_viper(
  mat = logFC_vector,
  network = regulons
)

# -----------------------------------------------
# Step 7E: Sort and Save TF Enrichment Results
# -----------------------------------------------

# Sort by absolute activity score
tf_results <- tf_activity[order(-abs(tf_activity$score)), ]

# Preview top TFs
head(tf_results, 10)

# Save to CSV
write.csv(tf_results,
          "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/tf_enrichment_results.csv",
          row.names = FALSE)


# -----------------------------------------------
# Step 7F: Barplot of Top 15 Transcription Factors (by activity score)
# -----------------------------------------------

library(ggplot2)

# Select top 15 by absolute activity score
top_tfs <- tf_results[order(-abs(tf_results$score)), ][1:15, ]

# Ensure correct TF order in the plot
top_tfs$source <- factor(top_tfs$source, levels = top_tfs$source[order(top_tfs$score)])

# Create the barplot
tf_barplot <- ggplot(top_tfs, aes(x = source, y = score, fill = score > 0)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +
  labs(title = "Top 15 Transcription Factors (VIPER TF Activity)",
       x = "Transcription Factor", y = "Activity Score (VIPER)") +
  scale_fill_manual(values = c("TRUE" = "#1B9E77", "FALSE" = "#D95F02")) +
  theme_minimal(base_size = 14)

# Display the plot
print(tf_barplot)

# Save to file
ggsave(filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/tf_activity_barplot.png",
       plot = tf_barplot,
       width = 10, height = 6, dpi = 300)




# -----------------------------------------------
# Step 8: Epigenetic Pathway Enrichment using GSEA
# -----------------------------------------------

# -----------------------------------------------
# Step 8A: Load Required Packages
# -----------------------------------------------


library(clusterProfiler)
library(msigdbr)
library(enrichplot)

# -----------------------------------------------
# Step 8B: Prepare Ranked Gene List (log2FC Vector)
# -----------------------------------------------
res_df_annotated <- readRDS("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/res_df_annotated.rds")

logFC_vector <- res_df_annotated$log2FoldChange
names(logFC_vector) <- res_df_annotated$SYMBOL
logFC_vector <- logFC_vector[is.finite(logFC_vector)]
logFC_vector <- sort(logFC_vector, decreasing = TRUE)

# -----------------------------------------------
# Step 8C: Get Epigenetic-Related Gene Sets from MSigDB (Reactome)
# -----------------------------------------------
msigdb_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

# Filter gene sets with epigenetic-related terms
epigenetic_sets <- msigdb_reactome[grep("CHROMATIN|METHYLATION|HISTONE|EPIGENETIC", msigdb_reactome$gs_name), ]

# Format to TERM2GENE
epigenetic_gene_sets <- split(epigenetic_sets$gene_symbol, epigenetic_sets$gs_name)
term2gene <- stack(epigenetic_gene_sets)[, c(2,1)]
colnames(term2gene) <- c("gene", "term")

# -----------------------------------------------
# Step 8D: Run GSEA with clusterProfiler
# -----------------------------------------------
gsea_epigenetic <- GSEA(
  geneList = logFC_vector,
  TERM2GENE = term2gene,
  pAdjustMethod = "BH",
  verbose = FALSE
)

# -----------------------------------------------
# Step 8E: Save GSEA Results to CSV
# -----------------------------------------------
write.csv(as.data.frame(gsea_epigenetic),
          "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/gsea_epigenetic_results.csv",
          row.names = FALSE)

# -----------------------------------------------
# Step 8F: Plot Top 10 Epigenetic Pathways (Dotplot)
# -----------------------------------------------

dotplot(gsea_epigenetic, showCategory = 10, title = "Epigenetic Pathway Enrichment (GSEA)")

# Save GSEA dotplot to PNG
ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/dotplot_gsea_epigenetic.png",
  plot = dotplot(gsea_epigenetic, showCategory = 10, title = "Epigenetic Pathway Enrichment (GSEA)"),
  width = 7, height = 4, dpi = 300
)

## Create barplot for top 10 enriched epigenetic pathways
# Convert GSEA results to data frame
gsea_df <- as.data.frame(gsea_epigenetic)

# Take top 10 based on absolute NES
top_gsea <- gsea_df %>%
  arrange(desc(abs(NES))) %>%
  head(10)

# Reorder terms by NES for plotting
top_gsea$Description <- factor(top_gsea$Description, levels = top_gsea$Description[order(top_gsea$NES)])

# Create barplot using ggplot2
barplot_gsea <- ggplot(top_gsea, aes(x = Description, y = NES, fill = p.adjust)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "#c11c84", high = "grey") +
  labs(
    title = "Epigenetic Pathway Enrichment (GSEA)",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "gray70", fill = NA)
  )

# Save the barplot as PNG
ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/barplot_gsea_epigenetic.png",
  plot = barplot_gsea,
  width = 10, height = 6, dpi = 300
)


# Ridgeplot
ridgeplot_gsea <- ridgeplot(gsea_epigenetic, showCategory = 10) +
  ggtitle("Epigenetic Pathway Enrichment (GSEA)")

# Save the plot
ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/ridgeplot_gsea_epigenetic.png",
  plot = ridgeplot_gsea,
  width = 10, height = 4, dpi = 300
)


# Ridgeplot with custom X-axis label
ridgeplot_gsea <- ridgeplot(gsea_epigenetic, showCategory = 10) +
  ggtitle("Epigenetic Pathway Enrichment (GSEA)") +
  xlab("Ranked Gene Position (Downregulated ←→ Upregulated)") +
  theme_minimal(base_size = 14)  # Optional: adjust text size

# Save the plot
ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/ridgeplot_gsea_epigenetic_labeled.png",
  plot = ridgeplot_gsea,
  width = 10, height = 6, dpi = 300
)




# MSigDB Reactome for Mitochondrial/Metabolic Gene Sets
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggplot2)

# -----------------------------------------------
# Step 8.2A: Load Required Packages
# -----------------------------------------------
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggplot2)
library(dplyr)

# -----------------------------------------------
# Step 8.2B: Prepare Ranked Gene List (log2FC Vector)
# -----------------------------------------------
res_df_annotated <- readRDS("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/res_df_annotated.rds")

logFC_vector <- res_df_annotated$log2FoldChange
names(logFC_vector) <- res_df_annotated$SYMBOL
logFC_vector <- logFC_vector[is.finite(logFC_vector)]
logFC_vector <- sort(logFC_vector, decreasing = TRUE)

# -----------------------------------------------
# Step 8.2C: Filter MSigDB Reactome for Mitochondrial/Metabolic Gene Sets
# -----------------------------------------------
msigdb_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

mito_sets <- msigdb_reactome[grep("MITOCHONDRIA|OXIDATIVE|RESPIRATION|ATP|TCA|CITRIC_ACID|ELECTRON_TRANSPORT|OXPHOS|METABOLISM", 
                                   msigdb_reactome$gs_name, ignore.case = TRUE), ]

# Format TERM2GENE
mito_gene_sets <- split(mito_sets$gene_symbol, mito_sets$gs_name)
term2gene_mito <- stack(mito_gene_sets)[, c(2, 1)]
colnames(term2gene_mito) <- c("gene", "term")

# -----------------------------------------------
# Step 8.2D: Run GSEA for Mitochondrial/Metabolism Pathways
# -----------------------------------------------
gsea_mito <- GSEA(
  geneList = logFC_vector,
  TERM2GENE = term2gene_mito,
  pAdjustMethod = "BH",
  verbose = FALSE
)

# -----------------------------------------------
# Step 8.2E: Save GSEA Results to CSV
# -----------------------------------------------
write.csv(as.data.frame(gsea_mito),
          "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/gsea_mitochondrial_results.csv",
          row.names = FALSE)

# -----------------------------------------------
# Step 8.2F: Plot Top 10 Mitochondrial Pathways (Dotplot)
# -----------------------------------------------
dotplot(gsea_mito, showCategory = 10, title = "Mitochondrial Pathway Enrichment (GSEA)")

ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/dotplot_gsea_mitochondria.png",
  plot = dotplot(gsea_mito, showCategory = 10, title = "Mitochondrial Pathway Enrichment (GSEA)"),
  width = 7, height = 4, dpi = 300
)

# -----------------------------------------------
# Step 8.2G: Barplot of Top 10 Mitochondrial Pathways
# -----------------------------------------------
gsea_mito_df <- as.data.frame(gsea_mito)

top_mito <- gsea_mito_df %>%
  arrange(desc(abs(NES))) %>%
  head(10)

top_mito$Description <- factor(top_mito$Description, levels = top_mito$Description[order(top_mito$NES)])

barplot_mito <- ggplot(top_mito, aes(x = Description, y = NES, fill = p.adjust)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "#4575b4", high = "grey") +
  labs(
    title = "Mitochondrial Pathway Enrichment (GSEA)",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "gray70", fill = NA)
  )

ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/barplot_gsea_mitochondria.png",
  plot = barplot_mito,
  width = 10, height = 6, dpi = 300
)

# -----------------------------------------------
# Step 8.2H: Ridgeplot of Top 10 Mitochondrial Pathways
# -----------------------------------------------
ridgeplot_mito <- ridgeplot(gsea_mito, showCategory = 10) +
  ggtitle("Mitochondrial Pathway Enrichment (GSEA)") +
  xlab("Ranked Gene Position (Downregulated ←→ Upregulated)") +
  theme_minimal(base_size = 14)

ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/ridgeplot_gsea_mitochondria.png",
  plot = ridgeplot_mito,
  width = 10, height = 6, dpi = 300
)





# -----------------------------------------------
# Step 9: Reactome Pathway Enrichment using GSEA (All Reactome Pathways)
# -----------------------------------------------

# -----------------------------------------------
# Step 9A: Load Required Packages
# -----------------------------------------------
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggplot2)
library(dplyr)

# -----------------------------------------------
# Step 9B: Prepare Ranked Gene List (log2FC Vector)
# -----------------------------------------------
res_df_annotated <- readRDS("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/res_df_annotated.rds")

logFC_vector <- res_df_annotated$log2FoldChange
names(logFC_vector) <- res_df_annotated$SYMBOL
logFC_vector <- logFC_vector[is.finite(logFC_vector)]
logFC_vector <- sort(logFC_vector, decreasing = TRUE)

# -----------------------------------------------
# Step 9C: Load All Reactome Gene Sets from MSigDB
# -----------------------------------------------
msigdb_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

# Format to TERM2GENE
reactome_gene_sets <- split(msigdb_reactome$gene_symbol, msigdb_reactome$gs_name)
term2gene <- stack(reactome_gene_sets)[, c(2, 1)]
colnames(term2gene) <- c("gene", "term")

# -----------------------------------------------
# Step 9D: Run GSEA with clusterProfiler
# -----------------------------------------------
gsea_reactome <- GSEA(
  geneList = logFC_vector,
  TERM2GENE = term2gene,
  pAdjustMethod = "BH",
  verbose = FALSE
)

# -----------------------------------------------
# Step 9E: Save GSEA Results to CSV
# -----------------------------------------------
write.csv(as.data.frame(gsea_reactome),
          "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/gsea_reactome_all_results.csv",
          row.names = FALSE)

# -----------------------------------------------
# Step 9F: Plot Top 15 Reactome Pathways (Dotplot)
# -----------------------------------------------
dotplot(gsea_reactome, showCategory = 15, title = "Reactome Pathway Enrichment (GSEA)")

# Save dotplot to PNG
ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/dotplot_gsea_reactome_all.png",
  plot = dotplot(gsea_reactome, showCategory = 15, title = "Reactome Pathway Enrichment (GSEA)"),
  width = 9, height = 10, dpi = 300
)

# -----------------------------------------------
# Step 9G: Custom Barplot for Top 15 Reactome Pathways
# -----------------------------------------------
gsea_df <- as.data.frame(gsea_reactome)

top_gsea <- gsea_df %>%
  arrange(desc(abs(NES))) %>%
  head(15)

top_gsea$Description <- factor(top_gsea$Description, levels = top_gsea$Description[order(top_gsea$NES)])

barplot_gsea <- ggplot(top_gsea, aes(x = Description, y = NES, fill = p.adjust)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "#c11c84", high = "grey") +
  labs(
    title = "Top Reactome Pathways (GSEA)",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "gray70", fill = NA)
  )

ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/barplot_gsea_reactome_all.png",
  plot = barplot_gsea,
  width = 10, height = 6, dpi = 300
)

# -----------------------------------------------
# Step 9H: Ridgeplot for Top 15 Reactome Pathways
# -----------------------------------------------
ridgeplot_gsea <- ridgeplot(gsea_reactome, showCategory = 15) +
  ggtitle("Reactome Pathway Enrichment (GSEA)")

ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/ridgeplot_gsea_reactome_all.png",
  plot = ridgeplot_gsea,
  width = 10, height = 12, dpi = 300
)

# Ridgeplot with labeled X-axis
ridgeplot_gsea <- ridgeplot(gsea_reactome, showCategory = 15) +
  ggtitle("Reactome Pathway Enrichment (GSEA)") +
  xlab("Ranked Gene Position (Downregulated ←→ Upregulated)") +
  theme_minimal(base_size = 14)

ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/ridgeplot_gsea_reactome_all_labeled.png",
  plot = ridgeplot_gsea,
  width = 10, height = 12, dpi = 300
)



# -----------------------------------------------
# Step 10: Hallmark Pathway Enrichment using GSEA (MSigDB H Collection)
# -----------------------------------------------

# Load required packages
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggplot2)
library(dplyr)

# Prepare ranked gene list (already done earlier, repeating here for clarity)
logFC_vector <- res_df_annotated$log2FoldChange
names(logFC_vector) <- res_df_annotated$SYMBOL
logFC_vector <- logFC_vector[is.finite(logFC_vector)]
logFC_vector <- sort(logFC_vector, decreasing = TRUE)

# Load Hallmark gene sets
msigdb_hallmark <- msigdbr(species = "Homo sapiens", category = "H")

# Format to TERM2GENE
hallmark_gene_sets <- split(msigdb_hallmark$gene_symbol, msigdb_hallmark$gs_name)
term2gene <- stack(hallmark_gene_sets)[, c(2, 1)]
colnames(term2gene) <- c("gene", "term")

# Run GSEA
gsea_hallmark <- GSEA(
  geneList = logFC_vector,
  TERM2GENE = term2gene,
  pAdjustMethod = "BH",
  verbose = FALSE
)

# Save GSEA results to CSV
write.csv(as.data.frame(gsea_hallmark),
          "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/gsea_hallmark_results.csv",
          row.names = FALSE)

# Dotplot
dotplot(gsea_hallmark, showCategory = 15, title = "Hallmark Pathway Enrichment (GSEA)")

ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/dotplot_gsea_hallmark.png",
  plot = dotplot(gsea_hallmark, showCategory = 10, title = "Hallmark Pathway Enrichment (GSEA)"),
  width = 9, height = 7, dpi = 300
)

# Custom barplot
gsea_df <- as.data.frame(gsea_hallmark)
top_gsea <- gsea_df %>% arrange(desc(abs(NES))) %>% head(15)
top_gsea$Description <- factor(top_gsea$Description, levels = top_gsea$Description[order(top_gsea$NES)])

barplot_gsea <- ggplot(top_gsea, aes(x = Description, y = NES, fill = p.adjust)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "#c11c84", high = "grey") +
  labs(
    title = "Top Hallmark Pathways (GSEA)",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "gray70", fill = NA)
  )

ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/barplot_gsea_hallmark.png",
  plot = barplot_gsea,
  width = 10, height = 6, dpi = 300
)

# Ridgeplot
ridgeplot_hallmark <- ridgeplot(gsea_hallmark, showCategory = 15) +
  ggtitle("Hallmark Pathway Enrichment (GSEA)") +
  xlab("Ranked Gene Position (Downregulated ←→ Upregulated)") +
  theme_minimal(base_size = 14)

ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/ridgeplot_gsea_hallmark_labeled.png",
  plot = ridgeplot_hallmark,
  width = 10, height = 6, dpi = 300
)




# -----------------------------------------------
# Step 11: GO Biological Process Enrichment using GSEA
# -----------------------------------------------

# -----------------------------------------------
# Step 11A: Load Required Packages
# -----------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# -----------------------------------------------
# Step 11B: Prepare Ranked Gene List (log2FC Vector)
# -----------------------------------------------
res_df_annotated <- readRDS("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/res_df_annotated.rds")

logFC_vector <- res_df_annotated$log2FoldChange
names(logFC_vector) <- res_df_annotated$SYMBOL
logFC_vector <- logFC_vector[is.finite(logFC_vector)]
logFC_vector <- sort(logFC_vector, decreasing = TRUE)

# -----------------------------------------------
# Step 11C: Run GSEA with GO:BP
# -----------------------------------------------
gsea_go <- gseGO(
  geneList = logFC_vector,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "SYMBOL",
  pAdjustMethod = "BH",
  verbose = FALSE
)

# -----------------------------------------------
# Step 11D: Save GSEA Results to CSV
# -----------------------------------------------
write.csv(as.data.frame(gsea_go),
          "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/gsea_go_results.csv",
          row.names = FALSE)

# -----------------------------------------------
# Step 11E: Dotplot for Top 15 GO Terms
# -----------------------------------------------
ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/dotplot_gsea_go.png",
  plot = dotplot(gsea_go, showCategory = 15, title = "GO:BP Enrichment (GSEA)"),
  width = 9, height = 10, dpi = 300
)

# -----------------------------------------------
# Step 11F: Barplot for Top 15 GO Terms
# -----------------------------------------------
gsea_df <- as.data.frame(gsea_go)
top_gsea <- gsea_df %>% arrange(desc(abs(NES))) %>% head(15)
top_gsea$Description <- factor(top_gsea$Description, levels = top_gsea$Description[order(top_gsea$NES)])

barplot_gsea <- ggplot(top_gsea, aes(x = Description, y = NES, fill = p.adjust)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "#c11c84", high = "grey") +
  labs(
    title = "Top GO:BP Pathways (GSEA)",
    x = "GO Term",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "gray70", fill = NA)
  )

ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/barplot_gsea_go.png",
  plot = barplot_gsea,
  width = 10, height = 6, dpi = 300
)

# -----------------------------------------------
# Step 11G: Ridgeplot with Label
# -----------------------------------------------
ridgeplot_go <- ridgeplot(gsea_go, showCategory = 15) +
  ggtitle("GO:BP Pathway Enrichment (GSEA)") +
  xlab("Ranked Gene Position (Downregulated ←→ Upregulated)") +
  theme_minimal(base_size = 14)

ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/ridgeplot_gsea_go.png",
  plot = ridgeplot_go,
  width = 10, height = 12, dpi = 300
)




# -----------------------------------------------
# Step 12: KEGG Pathway Enrichment using GSEA
# -----------------------------------------------

# -----------------------------------------------
# Step 12A: Load Required Packages
# -----------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# -----------------------------------------------
# Step 12B: Prepare Ranked Gene List (log2FC Vector)
# -----------------------------------------------
res_df_annotated <- readRDS("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/res_df_annotated.rds")

head(res_df_annotated)

# Map SYMBOL to ENTREZID
library(org.Hs.eg.db)
entrez_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = res_df_annotated$SYMBOL,
  keytype = "SYMBOL",
  columns = "ENTREZID"
)

# Merge Entrez IDs back into res_df_annotated
res_df_annotated <- merge(res_df_annotated, entrez_map, by = "SYMBOL")

# Now you can filter and continue
res_df_annotated <- res_df_annotated[!is.na(res_df_annotated$ENTREZID), ]

head(res_df_annotated)

logFC_vector <- res_df_annotated$log2FoldChange
names(logFC_vector) <- res_df_annotated$ENTREZID  # KEGG uses Entrez IDs
logFC_vector <- logFC_vector[is.finite(logFC_vector)]
logFC_vector <- sort(logFC_vector, decreasing = TRUE)

# -----------------------------------------------
# Step 12C: Run GSEA with KEGG
# -----------------------------------------------
gsea_kegg <- gseKEGG(
  geneList = logFC_vector,
  organism = "hsa",
  pAdjustMethod = "BH",
  verbose = FALSE
)

# -----------------------------------------------
# Step 12D: Save GSEA Results to CSV
# -----------------------------------------------
write.csv(as.data.frame(gsea_kegg),
          "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/gsea_kegg_results.csv",
          row.names = FALSE)

# -----------------------------------------------
# Step 12E: Dotplot for Top 15 KEGG Pathways
# -----------------------------------------------
ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/dotplot_gsea_kegg.png",
  plot = dotplot(gsea_kegg, showCategory = 15, title = "KEGG Pathway Enrichment (GSEA)"),
  width = 9, height = 8, dpi = 300
)

# -----------------------------------------------
# Step 12F: Barplot for Top 15 KEGG Pathways
# -----------------------------------------------
gsea_df <- as.data.frame(gsea_kegg)
top_gsea <- gsea_df %>% arrange(desc(abs(NES))) %>% head(15)
top_gsea$Description <- factor(top_gsea$Description, levels = top_gsea$Description[order(top_gsea$NES)])

barplot_gsea <- ggplot(top_gsea, aes(x = Description, y = NES, fill = p.adjust)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "#c11c84", high = "grey") +
  labs(
    title = "Top KEGG Pathways (GSEA)",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "gray70", fill = NA)
  )

ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/barplot_gsea_kegg.png",
  plot = barplot_gsea,
  width = 10, height = 6, dpi = 300
)

# -----------------------------------------------
# Step 12G: Ridgeplot for KEGG
# -----------------------------------------------
ridgeplot_kegg <- ridgeplot(gsea_kegg, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment (GSEA)") +
  xlab("Ranked Gene Position (Downregulated ←→ Upregulated)") +
  theme_minimal(base_size = 14)

ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/ridgeplot_gsea_kegg.png",
  plot = ridgeplot_kegg,
  width = 10, height = 12, dpi = 300
)



# -----------------------------------------------
# Step 13: GO Epigenetic Pathway Enrichment using GSEA
# -----------------------------------------------

# -----------------------------------------------
# Step 13A: Load Required Packages
# -----------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# -----------------------------------------------
# Step 13B: Prepare Ranked Gene List (log2FC Vector)
# -----------------------------------------------
res_df_annotated <- readRDS("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/res_df_annotated.rds")

logFC_vector <- res_df_annotated$log2FoldChange
names(logFC_vector) <- res_df_annotated$SYMBOL
logFC_vector <- logFC_vector[is.finite(logFC_vector)]
logFC_vector <- sort(logFC_vector, decreasing = TRUE)

# -----------------------------------------------
# Step 13C: Run GSEA with GO:BP
# -----------------------------------------------
gsea_go <- gseGO(
  geneList = logFC_vector,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "SYMBOL",
  pAdjustMethod = "BH",
  verbose = FALSE
)

# -----------------------------------------------
# Step 13D: Filter for Epigenetic-Related GO Terms
# -----------------------------------------------
gsea_go_df <- as.data.frame(gsea_go)
gsea_go_df <- gsea_go_df[grep("CHROMATIN|METHYLATION|HISTONE|EPIGENETIC", 
                              gsea_go_df$Description, ignore.case = TRUE), ]

# Save filtered result
write.csv(gsea_go_df,
          "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/gsea_go_epigenetic_filtered.csv",
          row.names = FALSE)

# -----------------------------------------------
# Step 13E: Plot Top Epigenetic GO Terms (Barplot)
# -----------------------------------------------
top_gsea <- gsea_go_df %>% arrange(desc(abs(NES))) %>% head(15)
top_gsea$Description <- factor(top_gsea$Description, levels = top_gsea$Description[order(top_gsea$NES)])

barplot_gsea <- ggplot(top_gsea, aes(x = Description, y = NES, fill = p.adjust)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "#c11c84", high = "grey") +
  labs(
    title = "Top Epigenetic GO:BP Pathways (GSEA)",
    x = "GO Term",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "gray70", fill = NA)
  )

ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/barplot_gsea_go_epigenetic.png",
  plot = barplot_gsea,
  width = 10, height = 6, dpi = 300
)

# -----------------------------------------------
# Step 13F: Ridgeplot for Epigenetic GO Terms
# -----------------------------------------------
ridgeplot_go <- ridgeplot(gsea_go_df, showCategory = 15) +
  ggtitle("Epigenetic GO:BP Pathway Enrichment (GSEA)") +
  xlab("Ranked Gene Position (Downregulated ←→ Upregulated)") +
  theme_minimal(base_size = 14)

ggsave(
  filename = "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/ridgeplot_gsea_go_epigenetic.png",
  plot = ridgeplot_go,
  width = 10, height = 6, dpi = 300
)
