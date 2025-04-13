# Gene Expression Analysis of Biguanide Response in Glioblastoma

**Project: GBM_biguanide_analysis**

## Overview

This repository contains an RNA-seq analysis of glioblastoma cell lines to investigate transcriptional differences between biguanide-sensitive and biguanide-resistant models. The study focuses on identifying differentially expressed genes (DEGs), enriched biological pathways, and transcription factor (TF) activities contributing to drug response, with a specific emphasis on mitochondrial metabolism and epigenetic regulation.

### Detailed Description

This project presents a comprehensive RNA-seq analysis of glioblastoma (GBM) cell lines to dissect the molecular mechanisms underlying differential responses to biguanide treatment. Biguanides (e.g., metformin, phenformin) are mitochondrial complex I inhibitors with potential anti-tumor effects, yet the transcriptional programs associated with resistance remain incompletely understood.

The analysis includes:

-   **Data preprocessing and normalization** using DESeq2\
-   **Differential expression analysis** to identify DEGs between sensitive and resistant lines\
-   **Over-representation analysis (ORA)** of up/downregulated genes using <GO:BP>, KEGG, and Reactome pathways\
-   **Gene Set Enrichment Analysis (GSEA)** with MSigDB Hallmark, Reactome, and mitochondrial/epigenetic terms\
-   **Transcription factor activity inference** using VIPER + DoRothEA regulons\
-   **Expression heatmaps** for top DEGs, epigenetic genes, and TF targets\
-   **Volcano plots** to visualize gene-level significance and fold change

The pipeline emphasizes the identification of: - Pathways involved in **mitochondrial metabolism and OXPHOS**\
- Genes regulating **chromatin, methylation, and histone modification**\
- Transcription factors with potential regulatory roles in biguanide resistance, such as **FOS**, **HNF4A**, and **ATF4**

All code, figures, and intermediate files are structured for reproducibility and downstream interpretation.

## Objectives

-   Identify DEGs between biguanide-sensitive and resistant glioblastoma lines\
-   Perform enrichment analyses: GO, KEGG, Reactome, Hallmark, and Epigenetic-related pathways\
-   Assess transcription factor activities using **DoRothEA** and **decoupleR**\
-   Visualize gene expression and pathway activity via heatmaps, dot plots, ridge plots, and volcano plots\
-   Highlight potential biomarkers and regulatory drivers of drug response

## Dataset

-   Raw counts: Paired-end RNA-seq counts generated using **featureCounts**
-   Samples:
    -   **Sensitive**: PDX22, PDX43, PDX59, HK281\
    -   **Resistant**: G83, PDX10

## Methodology

1.  **Preprocessing**
    -   Clean count matrix and metadata creation
    -   DESeq2 normalization and differential expression
2.  **QC and Visualization**
    -   PCA, sample distance heatmaps, MA plots, and dispersion estimates
3.  **DEG Identification**
    -   log2FC \> 1 or \< -1, FDR \< 0.05
4.  **Over-Representation Analysis (ORA)**
    -   **GO:BP**, **KEGG**, **Reactome** based on significant up/down genes
5.  **Gene Set Enrichment Analysis (GSEA)**
    -   **Hallmark**, **Reactome**, **GO:BP,** **KEGG**
    -   Epigenetic and mitochondrial pathway-focused GSEA using MSigDB
6.  **Transcription Factor (TF) Activity Analysis**
    -   VIPER algorithm via **decoupleR**
    -   Regulons from **DoRothEA** (A/B confidence only)
7.  **Expression Visualization**
    -   Heatmaps for top 50 DEGs, TF targets, epigenetic-related genes
    -   Volcano plots (DEGs, epigenetic, TFs, EnhancedVolcano-based)

## Key Findings

-   Enriched pathways in resistant lines include mitochondrial metabolism, chromatin regulation, and stress-response signaling\
-   TFs such as **HNF4A**, **FOS**, and **ATF4** show distinct regulatory roles between conditions\
-   Epigenetic regulators are significantly enriched in the top GSEA pathways and DE TFs

## Dependencies

-   **R 4.3.0+**
-   **DESeq2**
-   **clusterProfiler**, **enrichplot**, **ReactomePA**
-   **msigdbr**, **org.Hs.eg.db**
-   **decoupleR**, **dorothea**, **viper**
-   **pheatmap**, **ggplot2**, **ggrepel**, **EnhancedVolcano**

## Output Examples

-   **Top 50 DEGs**\
    `heatmap_top50_DEGs_02.png`, `volcano_top50_DEGs.png`

-   **TF Activity & Targets**\
    `tf_activity_barplot.png`, `tf_target_gene_expression_heatmap_grouped.png`

-   **Epigenetic Signatures**\
    `dotplot_gsea_epigenetic.png`, `heatmap_epigenetic_TFs.png`

-   **Mitochondrial Pathways**\
    `barplot_gsea_mitochondria.png`, `ridgeplot_gsea_mitochondria.png`

## Citation

If you use this code or analysis in your work, please cite this repository and the tools used:\
**DESeq2**, **clusterProfiler**, **ReactomePA**, **DoRothEA**, **decoupleR**, **EnhancedVolcano**, and **msigdbr**.

![](https://github.com/chingyaousf/Gene-Expression-Analysis-of-Biguanide-Response-in-Glioblastoma-GBM_biguanide_analysis-/blob/main/analysis/qc_plots/MA_plot_big_labels.png?raw=true){width="523"}

![](https://github.com/chingyaousf/Gene-Expression-Analysis-of-Biguanide-Response-in-Glioblastoma-GBM_biguanide_analysis-/blob/main/analysis/qc_plots/PCA_plot.png?raw=true){width="611"}

![](https://github.com/chingyaousf/Gene-Expression-Analysis-of-Biguanide-Response-in-Glioblastoma-GBM_biguanide_analysis-/blob/main/analysis/qc_plots/sample_distance_heatmap.png?raw=true){width="523"}
