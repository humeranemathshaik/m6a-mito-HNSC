

#result graphs and and full pipeline would be uploaded once this article gets published  




setwd('c:\\mito3')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(maftools)
library(pheatmap)

# ----------------------------------------------------------------------
#  Download & Prepare Transcriptome Data (Store in explicit variables)
# ----------------------------------------------------------------------
query_tcga <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query_tcga)

tcga_se <- GDCprepare(query_tcga)

# ----------------------------------------------------------------------
#  Extract Assay Data (tpm_unstrand or unstranded)
# ----------------------------------------------------------------------
# Check assay names to ensure index 4 or named assay is selected safely
# standard assay for TPM: assay(tcga_se, "tpm_unstrand")
tcga_tpm <- assay(tcga_se,"tpm_unstrand")

# Map Ensembl IDs to standard Gene Symbols
rownames(tcga_tpm) <- rowData(tcga_se)$gene_name

# Filter for Primary Tumor samples (barcode position 14-15 == "01")
tumor_idx <- which(substr(colnames(tcga_tpm), 14, 15) == "01")
tumor_tpm <- tcga_tpm[, tumor_idx]

# Remove duplicated gene symbols if present
tumor_tpm <- tumor_tpm[!duplicated(rownames(tumor_tpm)), ]

# ----------------------------------------------------------------------
#  Clean Gene Vectors (Sanitized HGNC Official Symbols)
# ----------------------------------------------------------------------
m6a_genes <- c(
  "METTL3", "METTL14", "WTAP", "KIAA1429", "RBM15", "RBM15B", 
  "ZC3H13", "METTL16", "CBLL1", "ZCCHC4", "METTL5", "FTO", 
  "ALKBH5", "ALKBH3", "YTHDF1", "YTHDF2", "YTHDF3", "YTHDC1", 
  "YTHDC2", "IGF2BP1", "IGF2BP2", "IGF2BP3", "HNRNPC", "HNRNPA2B1", 
  "EIF3A", "FMR1", "PRRC2A"
)

# Raw block parsing for mtrna_genes to properly split multiline string
mtrna_raw <- "PUS1\nQTRT1\nRCC1L\nREXO2\nRMND1\nRNASEH1\nRPUSD3\nRPUSD4\nSLIRP\nSUPV3L1\nTBRG4\nTEFM\nTFAM\nTFB1M\nTFB2M\nTHG1L\nTOP1MT\nTRIT1\nTRMT1\nTRMT10C\nTRMT2B\nTRMT5\nTRMT61B\nTRMU\nTRNT1\nTRUB2\nYBEY\nYRDC\nAARS2\nAURKAIP1\nC12orf65\nCARS2\nCHCHD1\nCOA3\nCOX14\nDAP3\nDARS2\nEARS2\nEXD2\nFARS2\nGADD45GIP1\nGARS1\nGATB\nGATC\nGFM1\nGFM2\nGTPBP10\nGUF1\nHARS2\nHEMK1\nIARS2\nKARS1\nLARS2\nMALSU1\nMARS2\nMETAP1D\nMIEF1\nMPV17L2\nMRPL1\nMRPL10\nMRPL11\nMRPL13\nMRPL14\nMRPL15\nMRPL16\nMRPL17\nMRPL18\nMRPL19\nMRPL2\nMRPL20\nMRPL21\nMRPL22\nMRPL23\nMRPL24\nMRPL27\nMRPL28\nMRPL3\nMRPL30\nMRPL32\nMRPL33\nMRPL34\nMRPL35\nMRPL36\nMRPL37\nMRPL38\nMRPL39\nMRPL4\nMRPL40\nMRPL41\nMRPL42\nMRPL43\nMRPL44\nMRPL45\nMRPL46\nMRPL48\nMRPL49\nMRPL50\nMRPL51\nMRPL52\nMRPL53\nMRPL54\nMRPL55\nMRPL57\nMRPL58\nMRPL9\nMRPS10\nMRPS11\nMRPS12\nMRPS14\nMRPS15\nMRPS16\nMRPS17\nMRPS18A\nMRPS18B\nMRPS18C\nMRPS2\nMRPS21\nMRPS22\nMRPS23\nMRPS24\nMRPS25\nMRPS26\nMRPS27\nMRPS28\nMRPS30\nMRPS31\nMRPS33\nMRPS34\nMRPS35\nMRPS36\nMRPS5\nMRPS6\nMRRF\nMTERF3\nMTERF4\nMTFMT\nMTG1\nMTG2\nMTIF2\nMTIF3\nMTRF1\nMTRF1L\nNARS2\nOXA1L\nPARS2\nPDF\nPTCD3\nPUSL1\nQRSL1\nRARS2\nRBFA\nSARS2\nTACO1\nTARS2\nTIMM21\nTSFM\nTUFM\nVARS2\nWARS2\nYARS2"
mtrna_genes <- unlist(strsplit(mtrna_raw, "\n"))

# ----------------------------------------------------------------------
#  Intersect with Expression Matrix Rownames
# ----------------------------------------------------------------------
valid_m6a <- intersect(m6a_genes, rownames(tumor_tpm))
valid_mtrna <- intersect(mtrna_genes, rownames(tumor_tpm))

message(paste("Found", length(valid_m6a), "out of", length(m6a_genes), "m6A genes."))
message(paste("Found", length(valid_mtrna), "out of", length(mtrna_genes), "mtRNA genes."))

# ----------------------------------------------------------------------
#  Calculate Pearson Correlation (|r| > 0.30)
# ----------------------------------------------------------------------
exp_mtrna <- log2(tumor_tpm[valid_mtrna, ] + 1)
exp_m6a   <- log2(tumor_tpm[valid_m6a, ] + 1)

cor_mat <- cor(t(exp_mtrna), t(exp_m6a), method = "pearson")

# Extract mitochondrial genes associated with at least one m6A regulator (|r| > 0.30)
correlated_mitogenes <- names(which(apply(cor_mat, 1, function(x) any(abs(x) > 0.30))))

message(paste("Identified", length(correlated_mitogenes), "mtRNA metabolic genes correlated with m6A regulators."))


library(pheatmap)

#  Filter the matrix to show only genes with at least one |r| > 0.30 correlation
sig_mtrna <- names(which(apply(cor_mat, 1, function(x) any(abs(x) > 0.30))))
plot_mat <- cor_mat[sig_mtrna, ]

#  Draw the heatmap
pheatmap(
  mat = plot_mat,
  color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100), # Blue-White-Red gradient
  breaks = seq(-1, 1, length.out = 101),                         # Force scale bounds from -1 to 1
  cluster_rows = TRUE,                                          # Cluster mtRNA genes
  cluster_cols = TRUE,                                          # Cluster m6A genes
  show_rownames = TRUE,                                         # Display mtRNA gene symbols
  show_colnames = TRUE,                                         # Display m6A gene symbols
  fontsize_row = 3,                                             # Adjust text size for long gene lists
  fontsize_col = 8,
  main = "Pearson Correlation: mtRNA Metabolism vs. m6A Regulators (TCGA-HNSC)"
)



#  Query Clinical Supplement data specifying the XML format
query_clin_TCGA <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Clinical",          
  data.type = "Clinical Supplement",
  data.format = "bcr xml"             
)


#  View the validated results summary table
output_query_TCGA <- getResults(query_clin_TCGA)
print(output_query_TCGA)


#  Download data - GDCdownload will now work smoothly
GDCdownload(query_clin_TCGA)


#  Prepare data - Load it into your R session
hnsc_clinical_tables <- GDCprepare(query_clin_TCGA)



#  Query  Copy Number Variation
query_cpn_TCGA <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number",
  sample.type = "Primary Tumor"
)

#  View the validated results summary table
output_query_TCGA <- getResults(query_cpn_TCGA)
print(query_cpn_TCGA)

GDCdownload(query_cpn_TCGA)




#  Prepare clinical data
hnsc_clinical_tables <- GDCprepare(query_clin_TCGA) 

#  Query Copy Number Variation (Added access parameter)
query_cpn_TCGA <- GDCquery(
  project = "TCGA-HNSC", 
  data.category = "Copy Number Variation", 
  data.type = "Gene Level Copy Number", 
  sample.type = "Primary Tumor",
  access = "open" # Required for gene-level copy number data
) 

# View and print the validated results summary table
output_query_TCGA <- getResults(query_cpn_TCGA) 
print(output_query_TCGA) # Corrected variable name from query_cpn_TCGA to output_query_TCGA

#  Download data
GDCdownload(query_cpn_TCGA)

#________Packgaes______________#
library("TCGAbiolinks") # bioconductor package
library("SummarizedExperiment") # bioconductor package
library("DESeq2") # bioconductor package
library("IHW") # bioconductor package
library("biomaRt") # bioconductor package
library("apeglm") # bioconductor package
 # CRAN package
library("RColorBrewer") # CRAN package
library("PCAtools") # bioconductor package
library(reshape2) # CRAN package

#  Install BiocManager if you don't have it already
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

#  Install DESeq2
BiocManager::install("DESeq2")
library("DESeq2")
BiocManager::install("PCAtools")
library("PCAtools")
BiocManager::install("apeglm")
library("apeglm")
BiocManager::install("biomaRt")
library("biomaRt")
BiocManager::install("IHW")
library("IHW")

dat <- GDCprepare(query = query_tcga, save = TRUE, save.filename = "exp.rda")

# exp matrix
rna <- as.data.frame(SummarizedExperiment::assay(dat))
# clinical data
clinical <- data.frame(dat@colData)

# Check how many tumor and control sample are there:
## data are stored under "definition" column in the clinical dataset.
table(clinical$definition)

# Also from sample id it is possible to count normal and tumor samples:
table(substr(colnames(rna),14,14))

#The count matrix (rna dataset) and the rows of the column data (clinical dataset) MUST be in the same order.

# to see whether all rows of clinical are present in rna datset
all(rownames(clinical) %in% colnames(rna))

# whether they are in the same order:
all(rownames(clinical) == colnames(rna))

#_______Making_Expression_Object__________#

#We will use the column “definition”, as the grouping variable for gene expression analysis. 

# replace spaces with "_" in levels of definition column
clinical$definition <-  gsub(" ", "_", clinical$definition)

# making the definition column as factor
clinical$definition <- as.factor(clinical$definition)
# relevling factor to ensure tumors would be compared to normal tissue.
levels(clinical$definition)
#
clinical$definition <- relevel(clinical$definition, ref = "Solid_Tissue_Normal")

# Making DESeqDataSet object which stores all experiment data
dds <- DESeqDataSetFromMatrix(countData = rna,
                              colData = clinical,
                              design = ~ definition)

# prefilteration: it is not necessary but recommended to filter out low expressed genes

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# data tranfromation
vsd <- vst(dds, blind=FALSE)

# making PC object
p <- pca(assay(vsd), metadata = colData(vsd), removeVar = 0.1)

# create PCA plot for PCA1 and PCA2
biplot(p, colby = "definition", lab = NULL, legendPosition = 'right', )


# Fol all top 10 possible combination 
pairsplot(p,
          components = getComponents(p, c(1:10)),
          triangle = TRUE, trianglelabSize = 12,
          hline = 0, vline = 0,
          pointSize = 0.4,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby = 'definition',
          title = 'Pairs plot', plotaxes = FALSE,
          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))


normal_idx <- substr(colnames(assay(vsd)),14,14) == "1"
n_sample <- assay(vsd)[, c(normal_idx) ]
colnames(n_sample) <- paste("NT_", substr(colnames(n_sample),1,12))

install.packages("RColorBrewer")


library(RColorBrewer)


sampleDists <- dist(t(n_sample))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(n_sample)
colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# Dissimilarity matrix calculation
sampleDists <- dist(t(n_sample))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(n_sample)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# heatmap visualization
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


# ----------------------------------------------------------------------
#  Differential Expression Analysis (DESeq2)
# ----------------------------------------------------------------------
dds <- DESeq(dds)

# Extract default contrast results (Tumor vs Solid_Tissue_Normal)
res <- results(dds, alpha = 0.05, lfcThreshold = 1.0)

# Convert results to dataframe and attach Gene Symbols directly from metadata
res_df <- as.data.frame(res)
res_df$ensembl_id <- rownames(res_df)

# Map Ensembl IDs to Symbols using rowData from original SummarizedExperiment
res_df$symbol <- rowData(tcga_se)$gene_name[match(res_df$ensembl_id, rowData(tcga_se)$gene_id)]

# ----------------------------------------------------------------------
#  Filter specifically for m6A & mtRNA Metabolism candidate genes
# ----------------------------------------------------------------------
candidate_genes <- c(valid_m6a, valid_mtrna)
res_candidates <- res_df %>% 
  filter(symbol %in% candidate_genes) %>% 
  arrange(padj, desc(abs(log2FoldChange)))

# Save all candidate DE results
write.csv(res_candidates, file = "m6A_mtRNA_DEG_results.csv", row.names = FALSE)

# ----------------------------------------------------------------------
#  Extract Top Dysregulated Target Genes for Plotting
# ----------------------------------------------------------------------
top_target_ids <- res_candidates$ensembl_id[1:min(25, nrow(res_candidates))]

# Extract VST normalized counts for top targets
dys_reg <- assay(vsd)[top_target_ids, ]
rownames(dys_reg) <- res_candidates$symbol[match(top_target_ids, res_candidates$ensembl_id)]

# Melt dataset for visualization
melted_norm_counts <- melt(as.matrix(dys_reg))
colnames(melted_norm_counts) <- c("gene", "samplename", "normalized_counts")

# Assign groups based on clinical metadata definition
normal_samples <- rownames(colData(vsd)[colData(vsd)$definition == "Solid_Tissue_Normal", ])
melted_norm_counts$group <- ifelse(melted_norm_counts$samplename %in% normal_samples, "Normal", "Tumor")

# Save normalized counts to CSV for downstream plotting
write.csv(melted_norm_counts, file = "top_candidate_counts_formatted.csv", row.names = FALSE)




# volcanoplot fix 
# Test against null hypothesis of 0 fold change
res <- results(dds, alpha = 0.05)

# Convert to dataframe and map symbols
res_df <- as.data.frame(res)
res_df$ensembl_id <- rownames(res_df)
res_df$symbol <- rowData(tcga_se)$gene_name[match(res_df$ensembl_id, rowData(tcga_se)$gene_id)]

# Filter for your candidate list
res_candidates <- res_df %>% 
  filter(symbol %in% c(valid_m6a, valid_mtrna))

library(ggplot2)
library(ggrepel)

# Define significance with log2FC > 0.58 (1.5-fold change) and padj < 0.05
res_candidates <- res_candidates %>%
  mutate(
    pathway = case_when(
      symbol %in% valid_m6a ~ "m6A Regulator",
      symbol %in% valid_mtrna ~ "mtRNA Metabolism"
    ),
    significance = case_when(
      padj < 0.05 & log2FoldChange > 0.58 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -0.58 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

# Plot with revised cutoffs
ggplot(res_candidates, aes(x = log2FoldChange, y = -log10(padj), color = significance, shape = pathway)) +
  geom_point(alpha = 0.8, size = 2.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8", "Not Significant" = "grey70")) +
  geom_text_repel(
    data = subset(res_candidates, padj < 0.05 & abs(log2FoldChange) > 0.58),
    aes(label = symbol),
    size = 3,
    max.overlaps = 15
  ) +
  theme_bw() +
  labs(
    title = "Volcano Plot: m6A & mtRNA Metabolism DEGs (TCGA-HNSC)",
    x = "log2(Fold Change)",
    y = "-log10(Adjusted P-Value)",
    color = "Expression",
    shape = "Gene Family"
  )




library(pheatmap)

# Extract significant DEG Ensembl IDs
sig_candidates <- res_candidates %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1.0)

# Subset VST matrix for significant candidate DEGs
heatmap_counts <- assay(vsd)[sig_candidates$ensembl_id, ]
rownames(heatmap_counts) <- sig_candidates$symbol

# Annotation for samples (Tumor vs Normal)
sample_annotation <- data.frame(
  Group = ifelse(colData(vsd)$definition == "Solid_Tissue_Normal", "Normal", "Tumor"),
  row.names = colnames(heatmap_counts)
)

# Annotation for gene classes (m6A vs mtRNA)
gene_annotation <- data.frame(
  Class = ifelse(rownames(heatmap_counts) %in% valid_m6a, "m6A Regulator", "mtRNA Metabolism"),
  row.names = rownames(heatmap_counts)
)

# Colors for annotations
ann_colors <- list(
  Group = c("Normal" = "#377EB8", "Tumor" = "#E41A1C"),
  Class = c("m6A Regulator" = "#4DAF4A", "mtRNA Metabolism" = "#984EA3")
)

# Draw Z-score heatmap
pheatmap(
  mat = heatmap_counts,
  scale = "row",
  annotation_col = sample_annotation,
  annotation_row = gene_annotation,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 6,
  color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
  main = "Expression Profile: Significant m6A & mtRNA DEGs"
)




# Plot top 12 significant candidate DEGs from melted_norm_counts
ggplot(melted_norm_counts, aes(x = group, y = normalized_counts, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Normal" = "#377EB8", "Tumor" = "#E41A1C")) +
  theme_classic() +
  labs(
    title = "Top Dysregulated Candidate Genes (Tumor vs Normal)",
    x = "Tissue Type",
    y = "VST Normalized Expression",
    fill = "Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold")
  )
