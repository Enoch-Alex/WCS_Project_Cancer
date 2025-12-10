# ==============================================================================
# Module 1: Differential Expression Analysis with DESeq2
# ==============================================================================
# Project: Comparing anti-androgen drug mechanisms in prostate cancer
# Goal: Identify differentially expressed genes for each drug treatment
# Method: DESeq2 with Wald tests and LFC shrinkage
# ==============================================================================

library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(VennDiagram)
library(UpSetR)

# Load processed data from Module 0
load("module0_deseq2_processed.RData")

cat("=" %>% rep(78) %>% paste(collapse = ""), "\n")
cat("MODULE 1: DIFFERENTIAL EXPRESSION ANALYSIS\n")
cat("=" %>% rep(78) %>% paste(collapse = ""), "\n\n")

# ==============================================================================
# Step 1: Define Contrasts of Interest
# ==============================================================================

cat("\nStep 1: Defining contrasts for analysis...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

groups <- levels(colData$Group)
cat("Available groups:\n")
print(groups)

contrasts_list <- list(
  "Bica_vs_Control_noDHT" = c("Group", "Bicalutamide_No_DHT", "Control_No_DHT"),
  "Enza_vs_Control_noDHT" = c("Group", "Enzalutamide_No_DHT", "Control_No_DHT"),
  "Apalu_vs_Control_noDHT" = c("Group", "Apalutamide_No_DHT", "Control_No_DHT"),
  "Bica_vs_Control_DHT" = c("Group", "Bicalutamide_Plus_DHT", "Control_Plus_DHT"),
  "Enza_vs_Control_DHT" = c("Group", "Enzalutamide_Plus_DHT", "Control_Plus_DHT"),
  "Apalu_vs_Control_DHT" = c("Group", "Apalutamide_Plus_DHT", "Control_Plus_DHT"),
  "DHT_effect_Control" = c("Group", "Control_Plus_DHT", "Control_No_DHT"),
  "DHT_effect_Bica" = c("Group", "Bicalutamide_Plus_DHT", "Bicalutamide_No_DHT"),
  "DHT_effect_Enza" = c("Group", "Enzalutamide_Plus_DHT", "Enzalutamide_No_DHT"),
  "DHT_effect_Apalu" = c("Group", "Apalutamide_Plus_DHT", "Apalutamide_No_DHT"),
  "Bica_vs_Enza_noDHT" = c("Group", "Bicalutamide_No_DHT", "Enzalutamide_No_DHT"),
  "Bica_vs_Apalu_noDHT" = c("Group", "Bicalutamide_No_DHT", "Apalutamide_No_DHT"),
  "Enza_vs_Apalu_noDHT" = c("Group", "Enzalutamide_No_DHT", "Apalutamide_No_DHT")
)

cat("\nTotal contrasts to analyze:", length(contrasts_list), "\n")

# ==============================================================================
# Step 2: Perform Differential Expression Testing
# ==============================================================================

perform_DE_analysis <- function(dds, contrast, contrast_name) {
  
  cat("\nAnalyzing:", contrast_name, "\n")
  
  # Get results
  res <- results(dds, 
                 contrast = contrast,
                 alpha = 0.05,
                 pAdjustMethod = "BH")
  
  # Apply LFC shrinkage
  res_shrink <- lfcShrink(dds,
                          contrast = contrast,
                          res = res,
                          type = "normal")
  
  # Convert to dataframe
  res_df <- as.data.frame(res_shrink) %>%
    rownames_to_column("GeneID") %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  # Add gene annotations
  if (exists("gene_annotations")) {
    # Convert gene_annotations GeneID to character to match res_df
    gene_annot_subset <- gene_annotations %>%
      select(GeneID, Symbol, Description) %>%
      mutate(GeneID = as.character(GeneID)) %>%  # Convert to character
      distinct()
    
    res_df <- res_df %>%
      left_join(gene_annot_subset, by = "GeneID")
  }
  
  # Create GeneName column (use Symbol if available, otherwise GeneID)
  if ("Symbol" %in% colnames(res_df)) {
    res_df$GeneName <- ifelse(is.na(res_df$Symbol) | res_df$Symbol == "", 
                              res_df$GeneID, 
                              res_df$Symbol)
  } else {
    res_df$GeneName <- res_df$GeneID
  }
  
  # Add regulation status
  res_df <- res_df %>%
    mutate(
      Regulation = case_when(
        is.na(padj) ~ "NS",
        padj < 0.05 & log2FoldChange > 1 ~ "UP",
        padj < 0.05 & log2FoldChange < -1 ~ "DOWN",
        TRUE ~ "NS"
      ),
      Contrast = contrast_name
    )
  
  # Summary statistics
  n_up <- sum(res_df$Regulation == "UP", na.rm = TRUE)
  n_down <- sum(res_df$Regulation == "DOWN", na.rm = TRUE)
  n_sig <- sum(res_df$padj < 0.05, na.rm = TRUE)
  
  cat("  Significant genes (padj < 0.05):", n_sig, "\n")
  cat("  Up-regulated (LFC > 1):", n_up, "\n")
  cat("  Down-regulated (LFC < -1):", n_down, "\n")
  
  summary_stats <- data.frame(
    Contrast = contrast_name,
    Total_Genes = nrow(res_df),
    Significant = n_sig,
    Upregulated = n_up,
    Downregulated = n_down,
    Mean_LFC_up = mean(res_df$log2FoldChange[res_df$Regulation == "UP"], na.rm = TRUE),
    Mean_LFC_down = mean(res_df$log2FoldChange[res_df$Regulation == "DOWN"], na.rm = TRUE)
  )
  
  return(list(
    results = res_df,
    summary = summary_stats,
    res_object = res_shrink
  ))
}

DE_results <- list()
DE_summaries <- list()

for (contrast_name in names(contrasts_list)) {
  result <- perform_DE_analysis(dds, contrasts_list[[contrast_name]], contrast_name)
  
  DE_results[[contrast_name]] <- result$results
  DE_summaries[[contrast_name]] <- result$summary
  
  write.csv(result$results,
            file = paste0("results/DE_", contrast_name, ".csv"),
            row.names = FALSE)
}

DE_summary_combined <- bind_rows(DE_summaries)
write.csv(DE_summary_combined, "results/DE_summary_all_contrasts.csv", row.names = FALSE)

cat("\n" , "=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("Differential expression analysis complete!\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
# ==============================================================================
# Step 3: Visualize Overall DE Results
# ==============================================================================

cat("\nStep 3: Creating summary visualizations...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

pdf("figures/08_DE_summary_barplot.pdf", width = 14, height = 8)

DE_summary_long <- DE_summary_combined %>%
  select(Contrast, Upregulated, Downregulated) %>%
  pivot_longer(cols = c(Upregulated, Downregulated),
               names_to = "Direction",
               values_to = "Count")

p1 <- ggplot(DE_summary_long, aes(x = Contrast, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  labs(title = "Differentially Expressed Genes Across All Contrasts",
       subtitle = "padj < 0.05, |log2FC| > 1",
       x = "Contrast",
       y = "Number of DE Genes",
       fill = "Regulation") +
  scale_fill_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8")) +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3)

print(p1)
dev.off()

cat("Summary barplot saved\n")

# ==============================================================================
# Step 4: Create Volcano Plots
# ==============================================================================

cat("\nStep 4: Creating volcano plots...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

create_volcano_plot <- function(DE_data, contrast_name, top_n = 15) {
  
  DE_data <- DE_data %>%
    filter(!is.na(padj), !is.na(log2FoldChange))
  
  top_genes <- DE_data %>%
    filter(Regulation != "NS") %>%
    arrange(padj) %>%
    head(top_n) %>%
    pull(GeneName)
  
  p <- ggplot(DE_data, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = Regulation), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("UP" = "#E41A1C", 
                                  "DOWN" = "#377EB8", 
                                  "NS" = "grey70")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey30") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
    theme_minimal() +
    labs(title = paste("Volcano Plot:", gsub("_", " ", contrast_name)),
         x = "log2 Fold Change",
         y = "-log10(adjusted p-value)",
         color = "Regulation") +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  
  if (length(top_genes) > 0) {
    label_data <- DE_data %>%
      filter(GeneName %in% top_genes)
    
    p <- p + geom_text_repel(data = label_data,
                             aes(label = GeneName),
                             size = 3,
                             max.overlaps = 20,
                             box.padding = 0.5,
                             point.padding = 0.3,
                             segment.color = "grey50")
  }
  
  return(p)
}

main_contrasts <- c("Bica_vs_Control_noDHT", 
                    "Enza_vs_Control_noDHT", 
                    "Apalu_vs_Control_noDHT",
                    "Bica_vs_Control_DHT",
                    "Enza_vs_Control_DHT",
                    "Apalu_vs_Control_DHT")

pdf("figures/09_volcano_plots_main.pdf", width = 16, height = 12)
par(mfrow = c(2, 3))

for (contrast_name in main_contrasts) {
  p <- create_volcano_plot(DE_results[[contrast_name]], contrast_name)
  print(p)
}

dev.off()

cat("Volcano plots saved\n")

# ==============================================================================
# Step 5: MA Plots
# ==============================================================================

cat("\nStep 5: Creating MA plots...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

pdf("figures/10_MA_plots.pdf", width = 16, height = 12)
par(mfrow = c(2, 3))

for (contrast_name in main_contrasts) {
  res_data <- DE_results[[contrast_name]]
  
  plot(res_data$baseMean, res_data$log2FoldChange,
       log = "x",
       pch = 20,
       cex = 0.5,
       col = ifelse(res_data$padj < 0.05, "red", "grey50"),
       xlab = "Mean of normalized counts",
       ylab = "log2 Fold Change",
       main = gsub("_", " ", contrast_name))
  abline(h = 0, col = "blue", lwd = 2)
  legend("topright", 
         legend = c("Significant", "Not significant"),
         col = c("red", "grey50"),
         pch = 20)
}

dev.off()

cat("MA plots saved\n")

# ==============================================================================
# Step 6: Compare Drug Responses
# ==============================================================================

cat("\nStep 6: Comparing gene sets across drugs...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

drug_contrasts <- c("Bica_vs_Control_noDHT", 
                    "Enza_vs_Control_noDHT", 
                    "Apalu_vs_Control_noDHT")

DE_genes_by_drug <- lapply(drug_contrasts, function(contrast) {
  DE_results[[contrast]] %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
    pull(GeneName)
})

names(DE_genes_by_drug) <- c("Bicalutamide", "Enzalutamide", "Apalutamide")

cat("Number of DE genes per drug:\n")
print(sapply(DE_genes_by_drug, length))

pdf("figures/11_venn_diagram_drugs.pdf", width = 10, height = 10)

venn.plot <- venn.diagram(
  x = DE_genes_by_drug,
  category.names = c("Bicalutamide", "Enzalutamide", "Apalutamide"),
  filename = NULL,
  output = TRUE,
  imagetype = "png",
  height = 2000,
  width = 2000,
  resolution = 300,
  lwd = 2,
  col = c("#E41A1C", "#377EB8", "#4DAF4A"),
  fill = c(alpha("#E41A1C", 0.3), alpha("#377EB8", 0.3), alpha("#4DAF4A", 0.3)),
  cex = 1.5,
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  main = "Differentially Expressed Genes by Drug",
  main.cex = 1.8
)

grid.draw(venn.plot)
dev.off()

cat("Venn diagram saved\n")

all_DE_genes <- unique(unlist(DE_genes_by_drug))

gene_presence <- data.frame(
  Gene = all_DE_genes,
  Bicalutamide = all_DE_genes %in% DE_genes_by_drug$Bicalutamide,
  Enzalutamide = all_DE_genes %in% DE_genes_by_drug$Enzalutamide,
  Apalutamide = all_DE_genes %in% DE_genes_by_drug$Apalutamide
)

gene_presence$N_Drugs <- rowSums(gene_presence[, -1])

gene_presence$Category <- case_when(
  gene_presence$N_Drugs == 3 ~ "All three drugs",
  gene_presence$N_Drugs == 2 & gene_presence$Bicalutamide & gene_presence$Enzalutamide ~ "Bica + Enza",
  gene_presence$N_Drugs == 2 & gene_presence$Bicalutamide & gene_presence$Apalutamide ~ "Bica + Apalu",
  gene_presence$N_Drugs == 2 & gene_presence$Enzalutamide & gene_presence$Apalutamide ~ "Enza + Apalu",
  gene_presence$N_Drugs == 1 & gene_presence$Bicalutamide ~ "Bicalutamide only",
  gene_presence$N_Drugs == 1 & gene_presence$Enzalutamide ~ "Enzalutamide only",
  gene_presence$N_Drugs == 1 & gene_presence$Apalutamide ~ "Apalutamide only"
)

write.csv(gene_presence, "results/gene_sharing_across_drugs.csv", row.names = FALSE)

cat("\nGene sharing summary:\n")
print(table(gene_presence$Category))

# ==============================================================================
# Step 7: UpSet Plot
# ==============================================================================

cat("\nStep 7: Creating UpSet plot...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

upset_data <- gene_presence %>%
  select(Gene, Bicalutamide, Enzalutamide, Apalutamide) %>%
  mutate(across(c(Bicalutamide, Enzalutamide, Apalutamide), as.numeric))

pdf("figures/12_upset_plot_drugs.pdf", width = 10, height = 6)
upset(upset_data,
      sets = c("Bicalutamide", "Enzalutamide", "Apalutamide"),
      order.by = "freq",
      main.bar.color = "#E41A1C",
      sets.bar.color = "#377EB8",
      text.scale = 1.5,
      point.size = 3.5,
      line.size = 1)
dev.off()

cat("UpSet plot saved\n")

# ==============================================================================
# Step 8: Heatmap of Top DE Genes
# ==============================================================================

cat("\nStep 8: Creating heatmap of top DE genes...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

top_genes_list <- lapply(drug_contrasts, function(contrast) {
  DE_results[[contrast]] %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    head(50) %>%
    pull(GeneID)
})

top_genes <- unique(unlist(top_genes_list))
cat("Total unique top genes:", length(top_genes), "\n")

heatmap_data <- vsd_mat[top_genes, ]
heatmap_data <- heatmap_data[complete.cases(heatmap_data), ]
heatmap_data_scaled <- t(scale(t(heatmap_data)))

if (exists("gene_annotations")) {
  gene_map <- gene_annotations %>%
    filter(GeneID %in% rownames(heatmap_data_scaled)) %>%
    select(GeneID, Symbol) %>%
    distinct()
  
  gene_names <- gene_map$Symbol
  names(gene_names) <- gene_map$GeneID
  
  new_rownames <- ifelse(rownames(heatmap_data_scaled) %in% names(gene_names),
                         gene_names[rownames(heatmap_data_scaled)],
                         rownames(heatmap_data_scaled))
  rownames(heatmap_data_scaled) <- new_rownames
}

annotation_col <- data.frame(
  Drug = colData$Drug,
  DHT = colData$DHT,
  row.names = colnames(heatmap_data_scaled)
)

ann_colors <- list(
  Drug = c(Control = "grey", 
           Bicalutamide = "#E41A1C", 
           Enzalutamide = "#377EB8", 
           Apalutamide = "#4DAF4A"),
  DHT = c(Plus_DHT = "#FF7F00", No_DHT = "#984EA3")
)

pdf("figures/13_heatmap_top_DE_genes.pdf", width = 12, height = 16)
pheatmap(heatmap_data_scaled,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 8,
         main = "Top Differentially Expressed Genes\n(Z-score normalized VST values)",
         border_color = NA)
dev.off()

cat("Heatmap saved\n")

# ==============================================================================
# Step 9: Direction of Change Analysis
# ==============================================================================

cat("\nStep 9: Analyzing direction of change across drugs...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

common_genes <- Reduce(intersect, DE_genes_by_drug)

cat("Genes common to all three drugs:", length(common_genes), "\n")

if (length(common_genes) > 0) {
  direction_comparison <- data.frame(Gene = common_genes)
  
  for (contrast_name in drug_contrasts) {
    drug_name <- gsub("_vs_Control_noDHT", "", contrast_name)
    
    fc_data <- DE_results[[contrast_name]] %>%
      filter(GeneName %in% common_genes) %>%
      select(GeneName, log2FoldChange, padj)
    
    colnames(fc_data) <- c("Gene", 
                           paste0("LFC_", drug_name), 
                           paste0("padj_", drug_name))
    
    direction_comparison <- direction_comparison %>%
      left_join(fc_data, by = "Gene")
  }
  
  direction_comparison <- direction_comparison %>%
    mutate(
      All_Same_Direction = sign(`LFC_Bica_vs_Control_noDHT`) == 
        sign(`LFC_Enza_vs_Control_noDHT`) &
        sign(`LFC_Enza_vs_Control_noDHT`) == 
        sign(`LFC_Apalu_vs_Control_noDHT`)
    )
  
  write.csv(direction_comparison, 
            "results/common_genes_direction_analysis.csv",
            row.names = FALSE)
  
  cat("Genes with same direction across all drugs:",
      sum(direction_comparison$All_Same_Direction, na.rm = TRUE), "\n")
  
  pdf("figures/14_fold_change_concordance.pdf", width = 12, height = 12)
  
  par(mfrow = c(2, 2))
  
  plot(direction_comparison$`LFC_Bica_vs_Control_noDHT`,
       direction_comparison$`LFC_Enza_vs_Control_noDHT`,
       pch = 20,
       col = alpha("blue", 0.6),
       xlab = "Bicalutamide log2FC",
       ylab = "Enzalutamide log2FC",
       main = "Bicalutamide vs Enzalutamide")
  abline(a = 0, b = 1, col = "red", lty = 2)
  abline(h = 0, v = 0, col = "grey", lty = 2)
  
  plot(direction_comparison$`LFC_Bica_vs_Control_noDHT`,
       direction_comparison$`LFC_Apalu_vs_Control_noDHT`,
       pch = 20,
       col = alpha("blue", 0.6),
       xlab = "Bicalutamide log2FC",
       ylab = "Apalutamide log2FC",
       main = "Bicalutamide vs Apalutamide")
  abline(a = 0, b = 1, col = "red", lty = 2)
  abline(h = 0, v = 0, col = "grey", lty = 2)
  
  plot(direction_comparison$`LFC_Enza_vs_Control_noDHT`,
       direction_comparison$`LFC_Apalu_vs_Control_noDHT`,
       pch = 20,
       col = alpha("blue", 0.6),
       xlab = "Enzalutamide log2FC",
       ylab = "Apalutamide log2FC",
       main = "Enzalutamide vs Apalutamide")
  abline(a = 0, b = 1, col = "red", lty = 2)
  abline(h = 0, v = 0, col = "grey", lty = 2)
  
  dev.off()
  
  cat("Fold change concordance plots saved\n")
}

# ==============================================================================
# Step 10: DHT Effect Analysis
# ==============================================================================

cat("\nStep 10: Analyzing DHT effects...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

dht_comparison <- data.frame(
  Drug = c("Control", "Bicalutamide", "Enzalutamide", "Apalutamide"),
  No_DHT = c(
    NA,
    sum(DE_results[["Bica_vs_Control_noDHT"]]$padj < 0.05, na.rm = TRUE),
    sum(DE_results[["Enza_vs_Control_noDHT"]]$padj < 0.05, na.rm = TRUE),
    sum(DE_results[["Apalu_vs_Control_noDHT"]]$padj < 0.05, na.rm = TRUE)
  ),
  Plus_DHT = c(
    NA,
    sum(DE_results[["Bica_vs_Control_DHT"]]$padj < 0.05, na.rm = TRUE),
    sum(DE_results[["Enza_vs_Control_DHT"]]$padj < 0.05, na.rm = TRUE),
    sum(DE_results[["Apalu_vs_Control_DHT"]]$padj < 0.05, na.rm = TRUE)
  ),
  DHT_Effect = c(
    sum(DE_results[["DHT_effect_Control"]]$padj < 0.05, na.rm = TRUE),
    sum(DE_results[["DHT_effect_Bica"]]$padj < 0.05, na.rm = TRUE),
    sum(DE_results[["DHT_effect_Enza"]]$padj < 0.05, na.rm = TRUE),
    sum(DE_results[["DHT_effect_Apalu"]]$padj < 0.05, na.rm = TRUE)
  )
)

write.csv(dht_comparison, "results/DHT_effect_summary.csv", row.names = FALSE)

cat("\nDHT effect summary:\n")
print(dht_comparison)

pdf("figures/15_DHT_effect_comparison.pdf", width = 10, height = 6)

dht_plot_data <- dht_comparison %>%
  pivot_longer(cols = c(No_DHT, Plus_DHT, DHT_Effect),
               names_to = "Condition",
               values_to = "DE_Genes") %>%
  filter(!is.na(DE_Genes))

ggplot(dht_plot_data, aes(x = Drug, y = DE_Genes, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Effect of DHT on Differential Expression",
       subtitle = "Number of significant DE genes (padj < 0.05)",
       x = "Treatment",
       y = "Number of DE Genes",
       fill = "Condition") +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

cat("DHT effect plot saved\n")

# ==============================================================================
# Step 11: Export Gene Lists
# ==============================================================================

cat("\nStep 11: Preparing gene lists for pathway analysis...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

dir.create("gene_lists", showWarnings = FALSE)

for (contrast_name in main_contrasts) {
  sig_genes <- DE_results[[contrast_name]] %>%
    filter(padj < 0.05) %>%
    select(GeneID, GeneName, log2FoldChange, padj, Regulation)
  
  write.csv(sig_genes,
            paste0("gene_lists/", contrast_name, "_significant.csv"),
            row.names = FALSE)
  
  up_genes <- sig_genes %>% filter(Regulation == "UP")
  write.csv(up_genes,
            paste0("gene_lists/", contrast_name, "_upregulated.csv"),
            row.names = FALSE)
  
  down_genes <- sig_genes %>% filter(Regulation == "DOWN")
  write.csv(down_genes,
            paste0("gene_lists/", contrast_name, "_downregulated.csv"),
            row.names = FALSE)
}

cat("Gene lists exported to gene_lists/ directory\n")

# ==============================================================================
# Step 12: Save All Results
# ==============================================================================

cat("\nStep 12: Saving all results...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

# Prepare list of objects to save (always needed)
objects_to_save <- c("DE_results", "DE_summary_combined", 
                     "dds", "vsd", "vsd_mat", "sample_info", "colData")

# Add optional objects if they exist
if (exists("gene_presence")) {
  objects_to_save <- c(objects_to_save, "gene_presence")
}

if (exists("direction_comparison")) {
  objects_to_save <- c(objects_to_save, "direction_comparison")
}

if (exists("dht_comparison")) {
  objects_to_save <- c(objects_to_save, "dht_comparison")
}

cat("Saving objects:", paste(objects_to_save, collapse = ", "), "\n")

# Save
save(list = objects_to_save,
     file = "module1_deseq2_results.RData")

cat("All results saved to: module1_deseq2_results.RData\n")

# ==============================================================================
# Final Summary
# ==============================================================================

summary_text <- paste0(
  "\n",
  "=" %>% rep(78) %>% paste(collapse = ""),
  "\n",
  "MODULE 1 COMPLETE - DIFFERENTIAL EXPRESSION ANALYSIS\n",
  "=" %>% rep(78) %>% paste(collapse = ""),
  "\n\n",
  "METHOD: DESeq2 with LFC shrinkage\n",
  "CONTRASTS ANALYZED: ", length(contrasts_list), "\n",
  "SIGNIFICANCE CUTOFF: padj < 0.05, |log2FC| > 1\n\n",
  "KEY FINDINGS:\n",
  "-" %>% rep(78) %>% paste(collapse = ""),
  "\n\n"
)

if (length(common_genes) > 0) {
  summary_text <- paste0(
    summary_text,
    "Drug-specific responses (no DHT):\n",
    "  Bicalutamide: ", length(DE_genes_by_drug$Bicalutamide), " genes\n",
    "  Enzalutamide: ", length(DE_genes_by_drug$Enzalutamide), " genes\n",
    "  Apalutamide: ", length(DE_genes_by_drug$Apalutamide), " genes\n\n",
    "Shared responses:\n",
    "  Common to all three drugs: ", length(common_genes), " genes\n\n"
  )
}

summary_text <- paste0(
  summary_text,
  "FILES CREATED:\n",
  "-" %>% rep(78) %>% paste(collapse = ""),
  "\n",
  "Results:\n",
  "  - DE_*.csv (individual contrast results)\n",
  "  - DE_summary_all_contrasts.csv\n",
  "  - gene_sharing_across_drugs.csv\n",
  "  - common_genes_direction_analysis.csv\n",
  "  - DHT_effect_summary.csv\n",
  "  - gene_lists/ (directory with gene lists)\n\n",
  "Figures:\n",
  "  - 08_DE_summary_barplot.pdf\n",
  "  - 09_volcano_plots_main.pdf\n",
  "  - 10_MA_plots.pdf\n",
  "  - 11_venn_diagram_drugs.pdf\n",
  "  - 12_upset_plot_drugs.pdf\n",
  "  - 13_heatmap_top_DE_genes.pdf\n",
  "  - 14_fold_change_concordance.pdf\n",
  "  - 15_DHT_effect_comparison.pdf\n\n",
  "NEXT STEP: Run Module 2 for pathway enrichment analysis\n",
  "=" %>% rep(78) %>% paste(collapse = ""),
  "\n"
)

cat(summary_text)
writeLines(summary_text, "results/MODULE1_SUMMARY.txt")

cat("\nModule 1 complete!\n\n")