# Project 3: 
# The intersection between adipose tissue depots and gender
# Author: Jakob Refsgaard Hartvig
# Date: 2025-01-31
# This program includes code and ideas developed with assistance from OpenAI's ChatGPT and Copilot Github. 
# Description: This script will analyze the intersection between adipose tissue depots and gender


# Load libraries, data and directory.
setwd("/Users/jakobhartvig/Desktop/Bio2_Eksamen")
library(readr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(DESeq2)
library(dplyr) 
library(tibble)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)


# Data loading
adipose_counts <- read.delim("adipose_tissue_counts.txt", sep = "\t", header = TRUE, row.names = 1) # Read the data and set the first column as row names
adipose_metadata <- read_delim("adipose_tissue_metadata.txt", delim = "\t", col_names = TRUE)


####################################################################
# Data distribution
####################################################################
# Generates a violin plot to visualize total RNA-seq counts across tissue depots and genders. 

# Combine total counts and metadata, add "GenderTissue" and "All" category in one step
total_counts_combined <- adipose_metadata %>%
  mutate(
    total_counts_adipose = colSums(adipose_counts)[Sample],
    GenderTissue = interaction(Gender, Tissue),
    All = "All"
  ) %>%
  bind_rows(mutate(., GenderTissue = All)) %>% 
  mutate(GenderTissue = factor(
    GenderTissue,
    levels = c("Female.Subcutaneous", "Female.Visceral",
               "Male.Subcutaneous", "Male.Visceral", "All")
  ))

# Includes boxplot, jittered points, and log-transformed y-axis for better distribution insight.
violin_total <- ggplot(total_counts_combined,
                       aes(x = GenderTissue, y = total_counts_adipose)) +
  geom_violin(trim = FALSE, aes(fill = GenderTissue)) +
  stat_summary(
    fun = median, geom = "point",
    color = "black", size = 2
  ) +
  labs(
    x = "Gender_Tissue",
    y = "Total Count"
  ) +
  theme_minimal() + 
  scale_y_log10() + 
  geom_violin(trim = FALSE, alpha = 0.1) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 12)
  )

# Print 
violin_total



####################################################################
# Differential gene expression, using DESeq2
####################################################################

# Convert Tissue and Gender columns in metadata to factors (ensures categorical treatment)
adipose_metadata <- adipose_metadata %>%
  mutate(
    Tissue = factor(Tissue, levels = c("Subcutaneous", "Visceral")),
    Gender = factor(Gender, levels = c("Female", "Male"))
  )

# Create DESeq2 dataset and run differential expression analysis
dds_adipose <- DESeqDataSetFromMatrix(
  countData = adipose_counts,  # Gene count matrix
  colData   = adipose_metadata,  # Metadata with Tissue and Gender
  design    = ~ Tissue + Gender  # Experimental design formula
)

dds_adipose <- DESeq(dds_adipose) # Performs differential expression analysis using the DESeq2 method. 

# Extract results for Tissue (Visceral vs. Subcutaneous) and Gender (Male vs. Female)
res_adipose_tissue <- results(dds_adipose, contrast = c("Tissue", "Visceral", "Subcutaneous"))
res_adipose_gender <- results(dds_adipose, contrast = c("Gender", "Male", "Female"))

# Extract significant DEGs padj < 0.05
degs_tissue <- res_adipose_tissue %>%
  as.data.frame() %>%
  filter(padj < 0.05)
head(degs_tissue)


# Extract significant DEGs padj < 0.05
degs_gender <- res_adipose_gender %>%
  as.data.frame() %>%
  filter(padj < 0.05)
head(degs_gender)


# Prepare results for visualization
# Format results for Male vs. Female and Visceral vs. Subcutaneous
# Adds significance classification and transformed p-values
# Gender: Male and female
res_gender <- as.data.frame(res_adipose_gender) %>%
  mutate(
    gene = rownames(res_adipose_gender),
    neg_log10_pvalue = -log10(pvalue),
    sig = case_when(
      padj < 0.05 & log2FoldChange >= 1  ~ "Upregulated",
      padj < 0.05 & log2FoldChange <= -1 ~ "Downregulated",
      TRUE                               ~ "NotSig"
    )
  )

# Adipose tissue: Visceral and Subcutaneous
res_tissue <- as.data.frame(res_adipose_tissue) %>%
  mutate(
    gene = rownames(res_adipose_tissue),
    neg_log10_pvalue = -log10(pvalue),
    sig = case_when(
      padj < 0.05 & log2FoldChange >= 1  ~ "Upregulated",
      padj < 0.05 & log2FoldChange <= -1 ~ "Downregulated",
      TRUE                               ~ "NotSig"
    )
  )


# Volcano plot function
# Generate a volcano plot
volcano_plot <- function(data, title) {
  ggplot(data, aes(x = log2FoldChange, y = neg_log10_pvalue)) +
    geom_point(aes(color = sig), alpha = 0.6) +  # Color points by significance
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "NotSig" = "grey")) +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 p-value") +
    theme_minimal() +
    theme(legend.title = element_blank(), legend.position = "top")
}


# Generate volcano plots

# Simple volcano plots
volcano_p1 <- volcano_plot(res_gender, title = "Volcano Plot: Male vs. Female")
volcano_p2 <- volcano_plot(res_tissue, title = "Volcano Plot: Visceral vs. Subcutaneous")

# Combine the two volcano plots
volcano_p1 + volcano_p2

# Add labels for top genes (e.g., significant and high fold changes)
p1_ggrepel <- volcano_plot(res_gender, title = "Volcano Plot: Male vs. Female") +
  geom_text_repel(
    data = res_gender %>% filter(padj < 0.05 & abs(log2FoldChange) > 1),
    aes(label = gene), max.overlaps = 10, size = 3
  )
p2_ggrepel <- volcano_plot(res_tissue, title = "Volcano Plot: Visceral vs. Subcutaneous") +
  geom_text_repel(
    data = res_tissue %>% filter(padj < 0.05 & abs(log2FoldChange) > 1),
    aes(label = gene), max.overlaps = 10, size = 3
  )

# Combine volcano plots with labels
p1_ggrepel + p2_ggrepel



####################################################################
# PCA
####################################################################
# Apply variance-stabilizing transformation to the count data
vsd_adipose <- vst(dds_adipose, blind = FALSE)

# PCA plot for tissue comparison (Subcutaneous vs. Visceral)
pca_1 <- plotPCA(vsd_adipose, intgroup = "Tissue") + 
  ggtitle("Adipose tissue, Sub. vs. Visc.")

# PCA plot for gender comparison (Female vs. Male)
pca_2 <- plotPCA(vsd_adipose, intgroup = "Gender") + 
  ggtitle("Adipose tissue, Female vs. Male")


# Add combined Group column to colData of dds_adipose directly
colData(dds_adipose)$Group <- with(colData(dds_adipose), interaction(Tissue, Gender, sep = "_"))
vsd_adipose <- vst(dds_adipose, blind = FALSE)

# PCA plot with combined Tissue and Gender groups
pca_3 <- plotPCA(vsd_adipose, intgroup = "Group") +
  ggtitle("PCA: Adipose tissue (Tissue and Gender Combined)")

# plotting of PCA
pca_1 / pca_2

pca_3

####################################################################
# Top 20 non-sig. genes
####################################################################
# Please note for the sake of space. TA = tissue adipose, GA = gender adipose

# Convert each results object to a data frame and retain key info
tissue_adipose <- as.data.frame(res_adipose_tissue) %>%
  rownames_to_column("Gene") %>%
  select(Gene, baseMean, log2FoldChange, padj) %>%
  rename(
    baseMean_TA = baseMean,
    lfc_TA      = log2FoldChange,
    padj_TA     = padj
  )

# Same for gender
gender_adipose <- as.data.frame(res_adipose_gender) %>%
  rownames_to_column("Gene") %>%
  select(Gene, baseMean, log2FoldChange, padj) %>%
  rename(
    baseMean_GA = baseMean,
    lfc_GA      = log2FoldChange,
    padj_GA     = padj
  )

# Join them by genes
combined_res_adipose <- tissue_adipose %>%
  left_join(gender_adipose, by = "Gene")

# Filter for non-significant genes
non_sig_genes_adipose <- combined_res_adipose %>%
  filter(
    !is.na(padj_TA),
    !is.na(padj_GA),
    padj_TA  >= 0.05,
    padj_GA >= 0.05
  )

# Calculate average baseMean across depot and gender
non_sig_genes_adipose <- non_sig_genes_adipose %>%
  mutate(
    baseMean_avg_adipose = (baseMean_TA + baseMean_GA) / 2
  )

# Select top 20 by average expression
top_20_non_sig_adipose <- non_sig_genes_adipose %>%
  arrange(desc(baseMean_avg_adipose)) %>%
  slice_head(n = 20)

# Print top non sig. genes
top_20_non_sig_adipose

####################
# Pathway Analysis
####################

# Extract top 20 non-significant genes for analysis
top_20_non_sig_genes <- top_20_non_sig_adipose$Gene

# Perform GO enrichment analysis (Biological Processes, p-value adjusted by BH method)
ego_adipose <- enrichGO(
  gene         = top_20_non_sig_genes,
  OrgDb        = org.Hs.eg.db, # data baase
  keyType      = "SYMBOL",   
  ont          = "BP",  # Biological process     
  pAdjustMethod= "BH", # padj
  pvalueCutoff = 0.05
)

# Subset to top 10 GO categories
ego_adipose_1_10 <- ego_adipose
ego_adipose_1_10@result <- ego_adipose@result[1:10, ]

# Dot plot for top 10 categories
nonsig_1_10 <- dotplot(ego_adipose_1_10, showCategory = 10) +
  ggtitle("GO Enrichment: Top 1-10, (Non-Sig. Genes)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 13),  
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
    legend.text = element_text(size = 12)
  )

# Subset to categories 11–20
ego_adipose_11_20 <- ego_adipose
ego_adipose_11_20@result <- ego_adipose@result[11:20, ]

# Dot plot for categories 11–20
nonsig_11_20 <- dotplot(ego_adipose_11_20, showCategory = 10) +
  ggtitle("GO Enrichment: Top 11–20 (Non-Sig. Genes)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 13),  
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
    legend.text = element_text(size = 12)
  )

# Combine the two dot plots
nonsig_1_10 + nonsig_11_20


####################################################################
# Heatmap for all differentially expressed genes 
####################################################################
# Z-score the VST-transformed counts

# Checking for NA values. 
sum(is.na(res_adipose_tissue$padj))
sum(is.na(res_adipose_gender$padj))


# Identify DE genes (removing rows with NA in padj)
de_genes <- unique(c(
  rownames(res_adipose_tissue[!is.na(res_adipose_tissue$padj) &
                                res_adipose_tissue$padj < 0.05, ]),
  rownames(res_adipose_gender[!is.na(res_adipose_gender$padj) &
                                res_adipose_gender$padj < 0.05, ])
))

# Extract VST counts for those DE genes
vsd_mat  <- assay(vsd_adipose)    # full VST-transformed matrix
vsd_de   <- vsd_mat[de_genes, ]   # subset to only DE genes

# z-scoring
vsd_de_z <- t(scale(t(vsd_de)))

# Column annotation
col_annotation <- HeatmapAnnotation(
  Tissue = factor(adipose_metadata$Tissue, 
                  levels = c("Subcutaneous", "Visceral")),
  Gender = factor(adipose_metadata$Gender, 
                  levels = c("Male", "Female")),
  col = list(
    Tissue = c("Subcutaneous" = "#84cb9f", "Visceral" = "#c0a5d8"),
    Gender  = c("Male" = "#56bcc2", "Female" = "#e77d72")
  ),
  annotation_name_side = "left"
)


# Create row clusters on the full matrix
row_dist <- dist(vsd_de_z)
row_hc   <- hclust(row_dist, method = "complete")
k        <- 4
row_clusters <- cutree(row_hc, k = k)
stopifnot(length(row_clusters) == nrow(vsd_de_z))

# Create a row annotation based on those 4 clusters
row_annotation <- rowAnnotation(
  Cluster = factor(row_clusters),
  col = list(
    Cluster = c(
      "1" = "#FF9999",
      "2" = "#9999FF",
      "3" = "#99FF99",
      "4" = "#FFFF99"
    )
  )
)

# Splitting the columns via the meta data
col_split <- data.frame(
  Tissue = factor(adipose_metadata$Tissue, 
                  levels = c("Subcutaneous", "Visceral")),
  Gender = factor(adipose_metadata$Gender, 
                  levels = c("Male", "Female"))
)


# Define a color scale for your Z‐scores (adjust the break points to your range)
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# plotting heatmap of z-scores, and the annotation made.
Heatmap(
  vsd_de_z,
  name            = "Z-score",
  col             = col_fun,
  top_annotation  = col_annotation,
  left_annotation = row_annotation,
  row_split       = row_clusters,
  column_split    = col_split,
  cluster_rows    = TRUE,
  cluster_columns = TRUE,
  show_row_names  = FALSE,
  show_column_names = FALSE
)



#################################################
# Pathway analysis per cluster
#################################################
go_results_list   <- list()

for (i in seq_len(k)) {
  # Genes in cluster i
  genes_in_cluster_i <- names(row_clusters)[row_clusters == i]
  
  # For KEGG enrichment, convert gene SYMBOLs to ENTREZIDs
  bitr_out <- bitr(
    genes_in_cluster_i, 
    fromType = "SYMBOL", 
    toType   = "ENTREZID", 
    OrgDb    = org.Hs.eg.db
  )
  
  # If fewer than 3 genes have valid ENTREZIDs, skip enrichment
  if (nrow(bitr_out) < 3) {
    message(sprintf("Skipping Cluster %d (only %d genes have ENTREZ IDs)", 
                    i, nrow(bitr_out)))
    next
  }
  
  # GO Enrichment (Biological Process)
  ego <- enrichGO(
    gene          = genes_in_cluster_i,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",   # If row names are SYMBOLs
    ont           = "BP",       # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )
  go_results_list[[paste0("Cluster_", i)]] <- ego
  
}


##################
# Inspect results
##################

# Plot the GO results for cluster 1 (top 10 terms)
go_1 <- dotplot(go_results_list[["Cluster_1"]], showCategory = 10) + 
  ggtitle("Cluster 1: GO Enrichment")

# Plot the GO results for cluster 2 (top 10 terms)
go_2 <- dotplot(go_results_list[["Cluster_2"]], showCategory = 10) + 
  ggtitle("Cluster 2: GO Enrichment")

# Plot the GO results for cluster 3 (top 10 terms)
go_3 <- dotplot(go_results_list[["Cluster_3"]], showCategory = 10) + 
  ggtitle("Cluster 3: GO Enrichment")

# Plot the GO results for cluster 4 (top 10 terms)
go_4 <- dotplot(go_results_list[["Cluster_4"]], showCategory = 10) + 
  ggtitle("Cluster 4: GO Enrichment")

# Combine the plots
go_1 + go_2 + go_3 + go_4



#################################################
# Boxplot of the 3 pathways 
#################################################
# From the GO enrichment results, the following pathways were chosen:
# rRNA metabolic process, cluster 1
# Positive regulation of cellular catabolic process, cluster 2
# Cytoplasmic translation, cluster 4


#############################
## CLUSTER 1: "rRNA metabolic process"
#############################

# Identify genes
res_1 <- go_results_list[["Cluster_1"]]@result
row_1 <- res_1[res_1$Description == "rRNA metabolic process", ]
genes_1_str <- row_1$geneID
genes_1_vec <- strsplit(genes_1_str, "/")[[1]]

# Subset expression
common_1 <- intersect(genes_1_vec, rownames(vsd_mat))
expr_1 <- vsd_mat[common_1, , drop = FALSE]

# Long format + metadata
df_1 <- as.data.frame(t(expr_1))
df_1$SampleID <- rownames(df_1)
df_1_long <- df_1 %>%
  gather("Gene", "Expression", -SampleID) %>%
  left_join(adipose_metadata, by=c("SampleID"="Sample"))

# Plot all genes
plot1_all <- ggplot(df_1_long, aes(x=Tissue, y=Expression, fill=Gender)) +
  geom_boxplot(position=position_dodge(width=0.8), outlier.shape=NA) +
  geom_jitter(position=position_dodge(width=0.8), alpha=0.5, size=1) +
  facet_wrap(~ Gene, scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Cluster_1: rRNA metabolic process - All Genes")

# Compute difference & plot top 20
df_1_means <- df_1_long %>%
  group_by(Gene, Tissue) %>%
  summarize(mean_expr = mean(Expression), .groups="drop") %>%
  pivot_wider(names_from=Tissue, values_from=mean_expr) %>%
  mutate(diff_in_means = Subcutaneous - Visceral) %>%
  arrange(desc(abs(diff_in_means)))

# top 20
top_n <- 20
top_1_genes <- df_1_means$Gene[1:min(top_n, nrow(df_1_means))]
df_1_top <- df_1_long %>% filter(Gene %in% top_1_genes)
df_1_top$Gene <- factor(df_1_top$Gene, levels=df_1_means$Gene)

# plotting
plot1_top <- ggplot(df_1_top, aes(x=Tissue, y=Expression, fill=Gender)) +
  geom_boxplot(position=position_dodge(width=0.8), outlier.shape=NA) +
  geom_jitter(position=position_dodge(width=0.8), alpha=0.5, size=1) +
  facet_wrap(~ Gene, scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Cluster_1: rRNA metabolic process - Top 20")

#############################
## CLUSTER 2: "positive regulation of cellular catabolic process"
#############################

# Identify genes
res_2 <- go_results_list[["Cluster_2"]]@result
row_2 <- res_2[res_2$Description == "positive regulation of cellular catabolic process", ]
genes_2_str <- row_2$geneID
genes_2_vec <- strsplit(genes_2_str, "/")[[1]]

go_results_list
row_2
# Subset expression
common_2 <- intersect(genes_2_vec, rownames(vsd_mat))
expr_2 <- vsd_mat[common_2, , drop=FALSE]

# Long format + metadata
df_2 <- as.data.frame(t(expr_2))
df_2$SampleID <- rownames(df_2)
df_2_long <- df_2 %>%
  gather("Gene", "Expression", -SampleID) %>%
  left_join(adipose_metadata, by=c("SampleID"="Sample"))

# Plot all genes
plot2_all <- ggplot(df_2_long, aes(x=Tissue, y=Expression, fill=Gender)) +
  geom_boxplot(position=position_dodge(width=0.8), outlier.shape=NA) +
  geom_jitter(position=position_dodge(width=0.8), alpha=0.5, size=1) +
  facet_wrap(~ Gene, scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Cluster_2: positive regulation of cellular catabolic process - All Genes")

# Compute difference & plot top 20
df_2_means <- df_2_long %>%
  group_by(Gene, Tissue) %>%
  summarize(mean_expr=mean(Expression), .groups="drop") %>%
  pivot_wider(names_from=Tissue, values_from=mean_expr) %>%
  mutate(diff_in_means = Subcutaneous - Visceral) %>%
  arrange(desc(abs(diff_in_means)))

# top 20 
top_2_genes <- df_2_means$Gene[1:min(top_n, nrow(df_2_means))]
df_2_top <- df_2_long %>% filter(Gene %in% top_2_genes)
df_2_top$Gene <- factor(df_2_top$Gene, levels=df_2_means$Gene)

# plotting
plot2_top <- ggplot(df_2_top, aes(x=Tissue, y=Expression, fill=Gender)) +
  geom_boxplot(position=position_dodge(width=0.8), outlier.shape=NA) +
  geom_jitter(position=position_dodge(width=0.8), alpha=0.5, size=1) +
  facet_wrap(~ Gene, scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Cluster_2: catabolic process - Top 20")

#############################
## CLUSTER 4: "cytoplasmic translation"
#############################

# Identify genes
res_4 <- go_results_list[["Cluster_4"]]@result
row_4 <- res_4[res_4$Description == "cytoplasmic translation", ]
genes_4_str <- row_4$geneID
genes_4_vec <- strsplit(genes_4_str, "/")[[1]]

# Subset expression
common_4 <- intersect(genes_4_vec, rownames(vsd_mat))
expr_4 <- vsd_mat[common_4, , drop=FALSE]

# Long format + metadata
df_4 <- as.data.frame(t(expr_4))
df_4$SampleID <- rownames(df_4)
df_4_long <- df_4 %>%
  gather("Gene", "Expression", -SampleID) %>%
  left_join(adipose_metadata, by=c("SampleID"="Sample"))

# Plot all genes
plot4_all <- ggplot(df_4_long, aes(x=Tissue, y=Expression, fill=Gender)) +
  geom_boxplot(position=position_dodge(width=0.8), outlier.shape=NA) +
  geom_jitter(position=position_dodge(width=0.8), alpha=0.5, size=1) +
  facet_wrap(~ Gene, scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Cluster_4: cytoplasmic translation - All Genes")

# Compute difference & plot top 20
df_4_means <- df_4_long %>%
  group_by(Gene, Tissue) %>%
  summarize(mean_expr=mean(Expression), .groups="drop") %>%
  pivot_wider(names_from=Tissue, values_from=mean_expr) %>%
  mutate(diff_in_means = Subcutaneous - Visceral) %>%
  arrange(desc(abs(diff_in_means)))

# top 20 cluster 4 genes
top_4_genes <- df_4_means$Gene[1:min(top_n, nrow(df_4_means))]
df_4_top <- df_4_long %>% filter(Gene %in% top_4_genes)
df_4_top$Gene <- factor(df_4_top$Gene, levels=df_4_means$Gene)

# plotting
plot4_top <- ggplot(df_4_top, aes(x=Tissue, y=Expression, fill=Gender)) +
  geom_boxplot(position=position_dodge(width=0.8), outlier.shape=NA) +
  geom_jitter(position=position_dodge(width=0.8), alpha=0.5, size=1) +
  facet_wrap(~ Gene, scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Cluster_4: cytoplasmic translation - Top 20 by")

#############################
## Done. 2 plots per cluster:
# plot1_all, plot1_top
# plot2_all, plot2_top
# plot4_all, plot4_top

print(plot1_all)
print(plot1_top)

print(plot2_all)
print(plot2_top)

print(plot4_all)
print(plot4_top)
#############################


