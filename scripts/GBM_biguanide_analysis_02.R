# -----------------------------------------------
# Visualizing Expression of Target Genes for Top Transcription Factors
# -----------------------------------------------

# -----------------------------------------------
# Step 0: Load Necessary Libraries
# -----------------------------------------------
library(dplyr)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)

# -----------------------------------------------
# Step 1: Identify Top TFs (already done)
# -----------------------------------------------

# Load VIPER TF enrichment results
tf_results <- read.csv("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/tf_enrichment_results.csv")

## Sort TFs by positive activity score (descending)
# Filter for TFs with positive activity scores
positive_tfs <- tf_results %>% 
  filter(score > 0) %>%
  arrange(-score)

# View top 10 positively enriched TFs
head(positive_tfs, 10)

# Select top 3 TFs with positive scores
top_tf_names <- positive_tfs$source[1:3]
top_tf_names

# ---------------------------------------------------------------------------
## Sort TFs by absolute activity score (descending)
tf_results <- tf_results[order(-abs(tf_results$score)), ]
# View top 10 TFs (for inspection)
head(tf_results, 10)

# Select top 3 TFs for visualization
top_tf_names <- tf_results$source[1:3]
top_tf_names


#------------------------------------------------------------------------------
## directly assign TF
top_tf_names <- c("HNF4A")
top_tf_names <- c("CEBPD")
top_tf_names <- c("REST")
top_tf_names <- c("ATF4")
top_tf_names <- c("FOS")

# -----------------------------------------------
# Step 2: Load and Prepare DoRothEA Regulons
# -----------------------------------------------

library(dorothea)

# Load human DoRothEA regulons and filter for high-confidence (A & B)
data(dorothea_hs, package = "dorothea")
regulons <- subset(dorothea_hs, confidence %in% c("A", "B"))

# Rename columns to match decoupleR format
colnames(regulons)[colnames(regulons) == "tf"] <- "source"
colnames(regulons)[colnames(regulons) == "target"] <- "target"
colnames(regulons)[colnames(regulons) == "mor"] <- "mor"
colnames(regulons)[colnames(regulons) == "confidence"] <- "likelihood"

# Preview regulon data
head(regulons)

# -----------------------------------------------
# Step 3: Extract Target Genes for Selected TFs
# -----------------------------------------------

# Filter DoRothEA regulons for the top TFs selected earlier
top_tf_targets <- regulons %>%
  filter(source %in% top_tf_names) %>%
  select(source, target) %>%
  distinct()


top_tf_targets <- regulons %>%
  dplyr::filter(source %in% top_tf_names) %>%
  dplyr::select(source, target) %>%
  dplyr::distinct()


# Preview: top TFs and their associated target genes
head(top_tf_targets)

# Optionally, list all unique target genes across the selected TFs
selected_target_genes <- unique(top_tf_targets$target)
length(selected_target_genes)  # how many target genes total


# -----------------------------------------------
# Step 4: Subset Expression Matrix for These Target Genes
# -----------------------------------------------

library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(pheatmap)

# Load normalized counts
expression_matrix <- read.csv("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/normalized_counts.csv",
                              row.names = 1)

# Clean Ensembl IDs (remove version numbers)
ensembl_ids_clean <- sub("\\..*", "", rownames(expression_matrix))
expression_matrix$ENSEMBL <- ensembl_ids_clean

# Map to gene symbols
symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ensembl_ids_clean,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
)

# Merge expression with symbol map
expression_matrix_annot <- merge(expression_matrix, symbol_map, by.x = "ENSEMBL", by.y = "ENSEMBL")

# Filter and clean
expression_matrix_annot <- expression_matrix_annot[!is.na(expression_matrix_annot$SYMBOL), ]
expression_matrix_annot <- expression_matrix_annot[!duplicated(expression_matrix_annot$SYMBOL), ]
rownames(expression_matrix_annot) <- expression_matrix_annot$SYMBOL

# Keep only expression columns
sample_cols <- c("PDX22", "PDX43", "PDX59", "G83", "PDX10", "HK281")
expression_matrix_clean <- expression_matrix_annot[, sample_cols]

# Subset to TF target genes
target_genes_present <- selected_target_genes[selected_target_genes %in% rownames(expression_matrix_clean)]
target_expr <- expression_matrix_clean[target_genes_present, ]

# Scale expression across genes
scaled_expr <- t(scale(t(target_expr)))
scaled_expr <- scaled_expr[apply(scaled_expr, 1, function(x) all(!is.na(x))), ]

# Create condition annotation (Resistant first)
sample_order <- c("G83", "PDX10", "PDX22", "PDX43", "PDX59", "HK281")
scaled_expr_ordered <- scaled_expr[, sample_order]
annotation_col <- data.frame(Condition = c("Resistant", "Resistant", "Sensitive", "Sensitive", "Sensitive", "Sensitive"))
rownames(annotation_col) <- sample_order

# Plot and save heatmap
pheatmap(scaled_expr_ordered,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 12,
         #fontsize_row = 15,
         #fontsize_col = 15,
         fontsize_row = 5,
         fontsize_col = 15,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Expression of TF(FOS) Target Genes",
         #main = "Expression of TF(ATF4) Target Genes",
         #main = "Expression of TF(HNF4A) Target Genes",
         #main = "Expression of TF(CEBPD) Target Genes",
         #main = "Expression of TF(REST) Target Genes",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

# Save plot to file
#png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/tf_target_gene_expression_heatmap_negative_CEPBD.png",
#png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/tf_target_gene_expression_heatmap_negative_REST.png",
#png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/tf_target_gene_expression_heatmap_negative_HNF4A.png",
#png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/tf_target_gene_expression_heatmap_negative_ATF4.png",
png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/tf_target_gene_expression_heatmap_negative_FOS.png",
    width = 1200, height = 1000, res = 150)
pheatmap(scaled_expr_ordered,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 12,
         #fontsize_row = 15,
         #fontsize_col = 15,
         fontsize_row = 5,
         fontsize_col = 15,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Expression of TF(FOS) Target Genes",
         #main = "Expression of TF(ATF4) Target Genes",
         #main = "Expression of TF(HNF4A) Target Genes",
         #main = "Expression of TF(CEBPD) Target Genes",
         #main = "Expression of TF(REST) Target Genes",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()






# Re-load the full merged matrix if needed
expression_matrix <- read.csv("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/normalized_counts.csv", 
                              row.names = 1)


head(expression_matrix)

ensembl_ids_clean <- sub("\\..*", "", rownames(expression_matrix))

head(ensembl_ids_clean)

library(AnnotationDbi)
library(org.Hs.eg.db)

symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ensembl_ids_clean,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
)

head(symbol_map)

expression_matrix$ENSEMBL <- ensembl_ids_clean

head(expression_matrix)

expression_matrix_annot <- merge(expression_matrix, symbol_map, by.x = "ENSEMBL", by.y = "ENSEMBL")

head(expression_matrix_annot)

# Check the merged result
dim(expression_matrix_annot)
head(expression_matrix_annot)



# -----------------------------------------------
# Finalize: Clean matrix using gene SYMBOLs as rownames
# -----------------------------------------------

# Step 1: Filter out rows with missing or duplicated SYMBOLs
expression_matrix_annot <- expression_matrix_annot[!is.na(expression_matrix_annot$SYMBOL), ]

head(expression_matrix_annot)

expression_matrix_annot <- expression_matrix_annot[!duplicated(expression_matrix_annot$SYMBOL), ]

head(expression_matrix_annot)

# Step 2: Set SYMBOLs as rownames
rownames(expression_matrix_annot) <- expression_matrix_annot$SYMBOL

head(expression_matrix_annot)

# Step 3: Keep only the expression columns (sample counts)
expression_matrix_clean <- expression_matrix_annot[, c("PDX22", "PDX43", "PDX59", "G83", "PDX10", "HK281")]

head(expression_matrix_clean)

# Step 4: Preview
dim(expression_matrix_clean)
head(rownames(expression_matrix_clean))


# Subset to the TF target genes that are present
target_genes_present <- selected_target_genes[selected_target_genes %in% rownames(expression_matrix_clean)]

head(target_genes_present)

# Extract their expression
target_expr <- expression_matrix_clean[target_genes_present, ]

head(target_expr)

# Final check
dim(target_expr)
head(target_expr)


# -----------------------------------------------
# Step 5: Scale Expression and Prepare Sample Annotations (Resistant vs Sensitive)
# -----------------------------------------------

# Z-score scale across genes (rows)
scaled_expr <- t(scale(t(target_expr)))

# Remove genes with all NAs (from constant expression)
scaled_expr <- scaled_expr[apply(scaled_expr, 1, function(x) all(!is.na(x))), ]

# Define correct condition order: match column order in scaled_expr
annotation_col <- data.frame(
  Condition = c("Sensitive", "Sensitive", "Sensitive", "Resistant", "Resistant", "Sensitive")
)
rownames(annotation_col) <- colnames(scaled_expr)

# Double check
annotation_col


# -----------------------------------------------
# Step 6: Plot Heatmap and Save
# -----------------------------------------------

library(pheatmap)

# Draw heatmap
pheatmap(scaled_expr,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 10,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Expression of TF Target Genes (Scaled)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

# Save to file
ggsave_path <- "/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/tf_target_gene_expression_heatmap.png"

png(ggsave_path, width = 1200, height = 1000, res = 150)
pheatmap(scaled_expr,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 10,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Expression of TF Target Genes (Scaled)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()


# -----------------------------------------------
# Step 6 (Updated): Plot Heatmap with Grouped Conditions
# -----------------------------------------------

library(pheatmap)

# Reorder columns: Resistant first, then Sensitive
sample_order <- c("G83", "PDX10", "PDX22", "PDX43", "PDX59", "HK281")
scaled_expr_ordered <- scaled_expr[, sample_order]

# Update annotation_col to match new order
annotation_col_ordered <- annotation_col[sample_order, , drop = FALSE]

# Plot and save
pheatmap(scaled_expr_ordered,
         annotation_col = annotation_col_ordered,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 5,
         fontsize_col = 10,
         cluster_rows = TRUE,
         cluster_cols = FALSE,  # <-- Disable column clustering
         main = "Expression of TF Target Genes (Grouped by Condition)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

# Save version
png("/home/cyang40/chingyao/GBM_biguanide_analysis/analysis/enrichment/plots/tf_target_gene_expression_heatmap_grouped.png",
    width = 1200, height = 1000, res = 150)
pheatmap(scaled_expr_ordered,
         annotation_col = annotation_col_ordered,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 5,
         fontsize_col = 10,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Expression of TF Target Genes (Grouped by Condition)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()
