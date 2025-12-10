# ==============================================================================
# Module 0: Data Loading and Preprocessing with DESeq2
# ==============================================================================
# Project: Comparing anti-androgen drug mechanisms in prostate cancer
# Dataset: GSE211781 - Human LNCaP prostate cancer cells
# Method: DESeq2 for differential expression analysis
# ==============================================================================

# Load required libraries
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(GEOquery)
library(ggrepel)

# Set working directory to where your data files are located
data_path <- "~/Library/CloudStorage/OneDrive-UniversityofCambridge/Document/PhD_Stuff/Applications/Computational_systems_biology_training/Comp_Sys_Bio2025-main/prostate_cancer_project_GSE211781"

setwd(data_path)

# Create output directories
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("pathway_analysis", showWarnings = FALSE)

cat("=" %>% rep(78) %>% paste(collapse = ""), "\n")
cat("MODULE 0: DATA LOADING WITH DESeq2\n")
cat("=" %>% rep(78) %>% paste(collapse = ""), "\n\n")

# ==============================================================================
# Step 1: Load Sample Metadata
# ==============================================================================

cat("\nStep 1: Loading sample metadata...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

gse <- getGEO(filename = "GSE211781_series_matrix.txt.gz", GSEMatrix = TRUE)
pData_raw <- pData(gse)

# Extract relevant information from the actual data structure
sample_info <- pData_raw %>%
  as.data.frame() %>%
  rownames_to_column("Sample_ID") %>%
  select(Sample_ID, title, description, contains("characteristics"), source_name_ch1) %>%
  mutate(
    # Extract treatment from characteristics_ch1.2
    Treatment_Raw = characteristics_ch1.2,
    Treatment_Clean = gsub("treatment: ", "", Treatment_Raw),
    
    # Parse drug information
    Drug = case_when(
      grepl("ARN|Apalutamide", Treatment_Clean, ignore.case = TRUE) ~ "Apalutamide",
      grepl("BIC|Bicalutamide", Treatment_Clean, ignore.case = TRUE) ~ "Bicalutamide",
      grepl("ENZ|Enzalutamide", Treatment_Clean, ignore.case = TRUE) ~ "Enzalutamide",
      grepl("VEH|Vehicle", Treatment_Clean, ignore.case = TRUE) ~ "Control",
      TRUE ~ "Unknown"
    ),
    
    # Parse DHT status
    DHT = case_when(
      grepl("DHT\\+", Treatment_Clean) ~ "Plus_DHT",
      grepl("VEH\\+", Treatment_Clean) ~ "No_DHT",
      TRUE ~ "No_DHT"
    ),
    
    # Extract replicate from title or description
    Replicate = case_when(
      grepl("rep 1|rep1|biol rep 1", title, ignore.case = TRUE) ~ "Rep1",
      grepl("rep 2|rep2|biol rep 2", title, ignore.case = TRUE) ~ "Rep2",
      grepl("rep 3|rep3|biol rep 3", title, ignore.case = TRUE) ~ "Rep3",
      TRUE ~ "Rep1"
    ),
    
    # Create grouping variables
    Group = paste(Drug, DHT, sep = "_"),
    Sample_Label = paste(Drug, DHT, Replicate, sep = "_")
  )

cat("\nSample Information (first 6 rows):\n")
print(sample_info %>% select(Sample_ID, Drug, DHT, Replicate, Group) %>% head())

cat("\nExperimental Design:\n")
design_table <- table(sample_info$Drug, sample_info$DHT)
print(design_table)

# Check if we have all expected groups
expected_drugs <- c("Control", "Bicalutamide", "Enzalutamide", "Apalutamide")
missing_drugs <- setdiff(expected_drugs, unique(sample_info$Drug))
if (length(missing_drugs) > 0) {
  cat("\nWARNING: Missing expected drugs:", paste(missing_drugs, collapse = ", "), "\n")
  cat("Available drugs:", paste(unique(sample_info$Drug), collapse = ", "), "\n")
}

write.csv(sample_info, "results/sample_metadata.csv", row.names = FALSE)


# ==============================================================================
# Step 2: Load Count Matrix
# ==============================================================================

cat("\nStep 2: Loading raw count data...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

counts_raw <- read.delim("GSE211781_raw_counts_GRCh38.p13_NCBI.tsv.gz",
                         header = TRUE,
                         row.names = 1,
                         check.names = FALSE)

cat("Count matrix: ", nrow(counts_raw), "genes x", ncol(counts_raw), "samples\n")
cat("Column names preview:\n")
print(head(colnames(counts_raw)))

# ==============================================================================
# Step 3: Load Gene Annotations
# ==============================================================================

cat("\nStep 3: Loading gene annotations...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

gene_annotations <- read.delim("Human.GRCh38.p13.annot.tsv.gz",
                               header = TRUE,
                               stringsAsFactors = FALSE,
                               comment.char = "#")

cat("Loaded annotations for", nrow(gene_annotations), "genes\n")
cat("Annotation columns:\n")
print(colnames(gene_annotations))

# ==============================================================================
# Step 4: Match Samples and Create DESeq2 Object
# ==============================================================================

cat("\nStep 4: Creating DESeqDataSet object...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

# Match sample order between count matrix and metadata
sample_info$GSM <- str_extract(sample_info$Sample_ID, "GSM[0-9]+")
count_colnames <- colnames(counts_raw)

cat("Count matrix column names:\n")
print(count_colnames)
cat("\nMetadata GSM IDs:\n")
print(sample_info$GSM)

# Check if column names contain GSM IDs
if (any(grepl("GSM", count_colnames))) {
  cat("\nMatching by GSM ID...\n")
  counts_gsm <- str_extract(count_colnames, "GSM[0-9]+")
  
  # Reorder sample_info to match counts_raw column order
  sample_info <- sample_info[match(counts_gsm, sample_info$GSM), ]
  
  # Check match
  cat("Sample matching successful:\n")
  match_check <- data.frame(
    Count_Col = count_colnames,
    GSM = counts_gsm,
    Metadata_GSM = sample_info$GSM,
    Drug = sample_info$Drug,
    DHT = sample_info$DHT
  )
  print(match_check)
  
} else {
  # If no GSM in column names, just use order
  cat("\nNo GSM in column names. Using sample order...\n")
  if (ncol(counts_raw) == nrow(sample_info)) {
    colnames(counts_raw) <- sample_info$Sample_ID
  } else {
    stop("ERROR: Number of samples in count matrix doesn't match metadata!")
  }
}

# Ensure counts are integers (required by DESeq2)
counts_matrix <- round(as.matrix(counts_raw))

# Verify no NAs in counts
if (any(is.na(counts_matrix))) {
  cat("WARNING: Found NA values in count matrix. Replacing with 0.\n")
  counts_matrix[is.na(counts_matrix)] <- 0
}

# Create sample data frame for DESeq2
colData <- data.frame(
  row.names = colnames(counts_matrix),
  Drug = factor(sample_info$Drug, levels = c("Control", "Bicalutamide", "Enzalutamide", "Apalutamide")),
  DHT = factor(sample_info$DHT, levels = c("No_DHT", "Plus_DHT")),
  Group = factor(sample_info$Group),
  Replicate = factor(sample_info$Replicate)
)

cat("\ncolData summary:\n")
print(colData)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = colData,
  design = ~ Group
)

cat("\n✓ DESeqDataSet created successfully\n")
cat("  Genes:", nrow(dds), "\n")
cat("  Samples:", ncol(dds), "\n")

# ==============================================================================
# Step 5: Pre-filtering
# ==============================================================================

cat("\nStep 5: Pre-filtering low-count genes...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

cat("Genes before filtering:", sum(keep) + sum(!keep), "\n")
cat("Genes after filtering:", sum(keep), "\n")
cat("Genes removed:", sum(!keep), "(", round(100*sum(!keep)/(sum(keep)+sum(!keep)), 1), "%)\n")

# ==============================================================================
# Step 6: Run DESeq2 Analysis
# ==============================================================================

cat("\nStep 6: Running DESeq2 analysis...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("This may take a few minutes...\n")

dds <- DESeq(dds)

cat("✓ DESeq2 analysis complete!\n")

# ==============================================================================
# Step 7: Quality Control - Size Factors
# ==============================================================================

cat("\nStep 7: Examining normalization (size factors)...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

size_factors <- sizeFactors(dds)
cat("\nSize factors:\n")
print(round(size_factors, 3))

pdf("figures/01_size_factors.pdf", width = 12, height = 6)
size_factor_df <- data.frame(
  Sample = names(size_factors),
  Size_Factor = size_factors,
  Group = colData$Group
)

ggplot(size_factor_df, aes(x = Sample, y = Size_Factor, fill = Group)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "DESeq2 Size Factors (Normalization)",
       x = "Sample",
       y = "Size Factor") +
  scale_fill_brewer(palette = "Set3")
dev.off()

cat("✓ Size factor plot saved\n")

# ==============================================================================
# Step 8: Quality Control - Dispersion Estimates
# ==============================================================================

cat("\nStep 8: Examining dispersion estimates...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

pdf("figures/02_dispersion_estimates.pdf", width = 10, height = 8)
plotDispEsts(dds, main = "Dispersion Estimates")
dev.off()

cat("✓ Dispersion plot saved\n")

# ==============================================================================
# Step 9: Variance Stabilizing Transformation
# ==============================================================================

cat("\nStep 9: Variance stabilizing transformation...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)
cat("VST complete. Matrix dimensions:", dim(vsd_mat)[1], "x", dim(vsd_mat)[2], "\n")

write.csv(vsd_mat, "results/normalized_vst_values.csv")
cat("✓ Normalized values saved\n")

# ==============================================================================
# Step 10: Principal Component Analysis
# ==============================================================================

cat("\nStep 10: Principal Component Analysis...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

pdf("figures/03_PCA_plot.pdf", width = 12, height = 10)

p1 <- plotPCA(vsd, intgroup = c("Drug", "DHT")) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "PCA of Samples (VST-transformed counts)") +
  scale_color_brewer(palette = "Dark2")

print(p1)

pca_full <- prcomp(t(vsd_mat))
pca_df <- as.data.frame(pca_full$x[, 1:5])
pca_df <- cbind(pca_df, colData)

variance_explained <- round(100 * summary(pca_full)$importance[2, 1:5], 1)

p2 <- ggplot(pca_df, aes(x = PC2, y = PC3, color = Drug, shape = DHT)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA: PC2 vs PC3",
       x = paste0("PC2 (", variance_explained[2], "%)"),
       y = paste0("PC3 (", variance_explained[3], "%)")) +
  scale_color_brewer(palette = "Dark2")

print(p2)

dev.off()

cat("✓ PCA plots saved\n")
write.csv(pca_df, "results/PCA_coordinates.csv", row.names = FALSE)

# ==============================================================================
# Step 11: Sample Distances and Clustering
# ==============================================================================

cat("\nStep 11: Sample distance matrix and clustering...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

sampleDists <- dist(t(vsd_mat))
sampleDistMatrix <- as.matrix(sampleDists)

annotation_col <- data.frame(
  Drug = colData$Drug,
  DHT = colData$DHT,
  row.names = rownames(colData)
)

ann_colors <- list(
  Drug = c(Control = "grey", 
           Bicalutamide = "#E41A1C", 
           Enzalutamide = "#377EB8", 
           Apalutamide = "#4DAF4A"),
  DHT = c(Plus_DHT = "#FF7F00", No_DHT = "#984EA3")
)

pdf("figures/04_sample_distances.pdf", width = 12, height = 10)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-Sample Distances",
         fontsize = 8)
dev.off()

cat("✓ Distance heatmap saved\n")

# ==============================================================================
# Step 12: Library Size Summary
# ==============================================================================

cat("\nStep 12: Library size and gene detection summary...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

qc_summary <- data.frame(
  Sample = colnames(dds),
  Library_Size = colSums(counts(dds)),
  Genes_Detected = colSums(counts(dds) > 0),
  Size_Factor = sizeFactors(dds),
  Group = colData$Group,
  Drug = colData$Drug,
  DHT = colData$DHT
)

write.csv(qc_summary, "results/QC_summary.csv", row.names = FALSE)

pdf("figures/05_library_sizes.pdf", width = 12, height = 6)
ggplot(qc_summary, aes(x = Sample, y = Library_Size/1e6, fill = Group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Library Sizes",
       x = "Sample",
       y = "Counts (millions)") +
  scale_fill_brewer(palette = "Set3")
dev.off()

cat("✓ QC plots saved\n")

# ==============================================================================
# Step 13: Save Processed Data
# ==============================================================================

cat("\nStep 13: Saving processed data...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

save(dds, vsd, vsd_mat, sample_info, colData, gene_annotations, pca_df,
     file = "module0_deseq2_processed.RData")

cat("✓ All data saved to: module0_deseq2_processed.RData\n")

# ==============================================================================
# Final Summary
# ==============================================================================

summary_text <- paste0(
  "\n",
  "=" %>% rep(78) %>% paste(collapse = ""),
  "\n",
  "MODULE 0 COMPLETE - DESeq2 ANALYSIS READY\n",
  "=" %>% rep(78) %>% paste(collapse = ""),
  "\n\n",
  "Dataset: GSE211781\n",
  "Method: DESeq2\n",
  "Total samples: ", ncol(dds), "\n",
  "Genes analyzed: ", nrow(dds), "\n",
  "Groups: ", length(unique(colData$Group)), "\n\n",
  "FILES CREATED:\n",
  "-" %>% rep(78) %>% paste(collapse = ""),
  "\n",
  "Results:\n",
  "  - sample_metadata.csv\n",
  "  - normalized_vst_values.csv\n",
  "  - PCA_coordinates.csv\n",
  "  - QC_summary.csv\n\n",
  "Figures:\n",
  "  - 01_size_factors.pdf\n",
  "  - 02_dispersion_estimates.pdf\n",
  "  - 03_PCA_plot.pdf\n",
  "  - 04_sample_distances.pdf\n",
  "  - 05_library_sizes.pdf\n\n",
  "NEXT STEP: Run Module 1 for differential expression analysis\n",
  "=" %>% rep(78) %>% paste(collapse = ""),
  "\n"
)

cat(summary_text)
writeLines(summary_text, "results/MODULE0_SUMMARY.txt")

