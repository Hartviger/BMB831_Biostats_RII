Project Title:
The Intersection Between Adipose Tissue Depots and Gender

Author:
Jakob Refsgaard Hartvig

Date:
2025-01-31

# BMB831_Biostats_RII
 Class project analyzing transcriptomic differences in adipose tissue depots (subcutaneous and visceral) and their variation by gender. Using bulk RNA-seq data from the GTEx Project, the code performs differential expression analysis, visualization, and pathway enrichment to explore gene expression patterns and their role in tissue function.

Key Steps:
Data Loading:
Import count and metadata files.
Exploratory Analysis:
Violin Plot: Visualize RNA-seq count distribution.
PCA: Assess clustering by tissue and gender.
Differential Expression:
Run DESeq2 to identify DEGs for tissue and gender comparisons (padj < 0.05).
Generate volcano plots for significant genes.
Pathway Analysis:
Perform GO enrichment on both significant and non-significant genes.
Visualize results with dot plots and boxplots.
Heatmap:
Cluster DE genes into four groups and display expression levels with a heatmap.


Dependencies:
readr, ggplot2, DESeq2, dplyr, ComplexHeatmap, clusterProfiler, org.Hs.eg.db
