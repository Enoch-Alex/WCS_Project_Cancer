# ==============================================================================
# Module 2: Pathway Enrichment Analysis
# ==============================================================================
# Project: Comparing anti-androgen drug mechanisms in prostate cancer
# Goal: Identify biological pathways affected by each drug
# Methods: Over-representation analysis (ORA) and Gene Set Enrichment Analysis (GSEA)
# ==============================================================================

library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(DOSE)
library(ReactomePA)
library(msigdbr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

# Load DE results from Module 1
load("module1_deseq2_results.RData")

cat("=" %>% rep(78) %>% paste(collapse = ""), "\n")
cat("MODULE 2: PATHWAY ENRICHMENT ANALYSIS\n")
cat("=" %>% rep(78) %>% paste(collapse = ""), "\n\n")

# Create output directory
dir.create("pathway_analysis", showWarnings = FALSE)

# ==============================================================================
# Step 1: Gene ID Conversion (Symbol to Entrez)
# ==============================================================================

cat("\nStep 1: Converting gene symbols to Entrez IDs...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

convert_to_entrez <- function(gene_symbols) {
  entrez_ids <- mapIds(org.Hs.eg.db,
                       keys = gene_symbols,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")
  
  conversion_df <- data.frame(
    Symbol = gene_symbols,
    EntrezID = as.character(entrez_ids),
    stringsAsFactors = FALSE
  )
  
  conversion_df <- conversion_df %>%
    filter(!is.na(EntrezID))
  
  success_rate <- nrow(conversion_df) / length(gene_symbols) * 100
  cat("  Converted:", nrow(conversion_df), "/", length(gene_symbols), 
      "genes (", round(success_rate, 1), "%)\n")
  
  return(conversion_df)
}

# ==============================================================================
# Step 2: Prepare Pathway Databases
# ==============================================================================

cat("\nStep 2: Loading pathway databases from MSigDB...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

# WikiPathways
wikipathways_human <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS")
wikipathways_list <- split(wikipathways_human$entrez_gene, wikipathways_human$gs_name)
cat("  WikiPathways:", length(wikipathways_list), "pathways\n")

# Reactome
reactome_human <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
reactome_list <- split(reactome_human$entrez_gene, reactome_human$gs_name)
cat("  Reactome:", length(reactome_list), "pathways\n")

# BIOCARTA (alternative to KEGG)
biocarta_human <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA")
biocarta_list <- split(biocarta_human$entrez_gene, biocarta_human$gs_name)
cat("  BIOCARTA:", length(biocarta_list), "pathways\n")

# PID (Pathway Interaction Database - another KEGG alternative)
pid_human <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:PID")
pid_list <- split(pid_human$entrez_gene, pid_human$gs_name)
cat("  PID:", length(pid_list), "pathways\n")

# Hallmark
hallmark_human <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_human$entrez_gene, hallmark_human$gs_name)
cat("  Hallmark:", length(hallmark_list), "pathways\n")

# Gene Ontology
go_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
go_bp_list <- split(go_bp$entrez_gene, go_bp$gs_name)
cat("  GO Biological Process:", length(go_bp_list), "terms\n")

go_mf <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF")
go_mf_list <- split(go_mf$entrez_gene, go_mf$gs_name)
cat("  GO Molecular Function:", length(go_mf_list), "terms\n")

go_cc <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC")
go_cc_list <- split(go_cc$entrez_gene, go_cc$gs_name)
cat("  GO Cellular Component:", length(go_cc_list), "terms\n")

# ==============================================================================
# Step 2: Prepare Pathway Databases
# ==============================================================================

cat("\nStep 2: Loading pathway databases...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

# WikiPathways
wikipathways_human <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS")
wikipathways_list <- split(wikipathways_human$entrez_gene, wikipathways_human$gs_name)
cat("  WikiPathways:", length(wikipathways_list), "pathways\n")

# Reactome
reactome_human <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
reactome_list <- split(reactome_human$entrez_gene, reactome_human$gs_name)
cat("  Reactome:", length(reactome_list), "pathways\n")

# KEGG - accessed directly via clusterProfiler (requires internet)
cat("  KEGG: Will be accessed directly via enrichKEGG() function\n")

# Hallmark
hallmark_human <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_human$entrez_gene, hallmark_human$gs_name)
cat("  Hallmark:", length(hallmark_list), "pathways\n")

# Gene Ontology
go_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
go_bp_list <- split(go_bp$entrez_gene, go_bp$gs_name)
cat("  GO Biological Process:", length(go_bp_list), "terms\n")

# ==============================================================================
# Step 3: Over-Representation Analysis (ORA)
# ==============================================================================

cat("\nStep 3: Running over-representation analysis...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

run_ORA_analysis <- function(DE_data, contrast_name, background_genes) {
  
  cat("\nAnalyzing:", contrast_name, "\n")
  
  # Get significant genes
  sig_genes <- DE_data %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
    pull(GeneName)
  
  if (length(sig_genes) < 5) {
    cat("  Too few significant genes (", length(sig_genes), "). Skipping.\n")
    return(NULL)
  }
  
  cat("  Significant genes:", length(sig_genes), "\n")
  
  # Convert to Entrez IDs
  sig_entrez <- convert_to_entrez(sig_genes)
  bg_entrez <- convert_to_entrez(background_genes)
  
  if (nrow(sig_entrez) < 5) {
    cat("  Too few genes after ID conversion. Skipping.\n")
    return(NULL)
  }
  
  results_list <- list()
  
  # WikiPathways
  tryCatch({
    wiki_enrich <- enricher(
      gene = sig_entrez$EntrezID,
      universe = bg_entrez$EntrezID,
      TERM2GENE = data.frame(
        term = rep(names(wikipathways_list), sapply(wikipathways_list, length)),
        gene = unlist(wikipathways_list)
      ),
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 5,
      maxGSSize = 500
    )
    
    if (!is.null(wiki_enrich) && nrow(wiki_enrich@result) > 0) {
      results_list[["WikiPathways"]] <- wiki_enrich@result %>%
        mutate(Database = "WikiPathways", Contrast = contrast_name)
      cat("  WikiPathways: ", nrow(wiki_enrich@result), "pathways enriched\n")
    }
  }, error = function(e) {
    cat("  WikiPathways: Error -", e$message, "\n")
  })
  
  # KEGG
  # KEGG (using enrichKEGG function instead)
  tryCatch({
    kegg_enrich <- enrichKEGG(
      gene = sig_entrez$EntrezID,
      organism = "hsa",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      universe = bg_entrez$EntrezID,
      minGSSize = 5,
      maxGSSize = 500
    )
    
    if (!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
      results_list[["KEGG"]] <- kegg_enrich@result %>%
        mutate(Database = "KEGG", Contrast = contrast_name)
      cat("  KEGG: ", nrow(kegg_enrich@result), "pathways enriched\n")
    }
  }, error = function(e) {
    cat("  KEGG: Error -", e$message, "\n")
    cat("  Note: KEGG requires internet connection\n")
  })
  
  # Reactome
  tryCatch({
    reactome_enrich <- enricher(
      gene = sig_entrez$EntrezID,
      universe = bg_entrez$EntrezID,
      TERM2GENE = data.frame(
        term = rep(names(reactome_list), sapply(reactome_list, length)),
        gene = unlist(reactome_list)
      ),
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 5,
      maxGSSize = 500
    )
    
    if (!is.null(reactome_enrich) && nrow(reactome_enrich@result) > 0) {
      results_list[["Reactome"]] <- reactome_enrich@result %>%
        mutate(Database = "Reactome", Contrast = contrast_name)
      cat("  Reactome: ", nrow(reactome_enrich@result), "pathways enriched\n")
    }
  }, error = function(e) {
    cat("  Reactome: Error -", e$message, "\n")
  })
  
  # Hallmark
  tryCatch({
    hallmark_enrich <- enricher(
      gene = sig_entrez$EntrezID,
      universe = bg_entrez$EntrezID,
      TERM2GENE = data.frame(
        term = rep(names(hallmark_list), sapply(hallmark_list, length)),
        gene = unlist(hallmark_list)
      ),
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 5,
      maxGSSize = 500
    )
    
    if (!is.null(hallmark_enrich) && nrow(hallmark_enrich@result) > 0) {
      results_list[["Hallmark"]] <- hallmark_enrich@result %>%
        mutate(Database = "Hallmark", Contrast = contrast_name)
      cat("  Hallmark: ", nrow(hallmark_enrich@result), "pathways enriched\n")
    }
  }, error = function(e) {
    cat("  Hallmark: Error -", e$message, "\n")
  })
  
  # GO Biological Process
  tryCatch({
    go_bp_enrich <- enricher(
      gene = sig_entrez$EntrezID,
      universe = bg_entrez$EntrezID,
      TERM2GENE = data.frame(
        term = rep(names(go_bp_list), sapply(go_bp_list, length)),
        gene = unlist(go_bp_list)
      ),
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 5,
      maxGSSize = 500
    )
    
    if (!is.null(go_bp_enrich) && nrow(go_bp_enrich@result) > 0) {
      results_list[["GO_BP"]] <- go_bp_enrich@result %>%
        mutate(Database = "GO_BP", Contrast = contrast_name)
      cat("  GO BP: ", nrow(go_bp_enrich@result), "terms enriched\n")
    }
  }, error = function(e) {
    cat("  GO BP: Error -", e$message, "\n")
  })
  
  return(results_list)
}

# Define background (all genes tested)
background_genes <- unique(unlist(lapply(DE_results, function(x) {
  x %>% filter(!is.na(padj)) %>% pull(GeneName)
})))

cat("Background universe:", length(background_genes), "genes\n")

# Run ORA for main contrasts
main_contrasts <- c("Bica_vs_Control_noDHT", 
                    "Enza_vs_Control_noDHT", 
                    "Apalu_vs_Control_noDHT",
                    "Bica_vs_Control_DHT",
                    "Enza_vs_Control_DHT",
                    "Apalu_vs_Control_DHT")

all_enrichments <- list()

for (contrast_name in main_contrasts) {
  enrichment_results <- run_ORA_analysis(
    DE_results[[contrast_name]], 
    contrast_name, 
    background_genes
  )
  
  if (!is.null(enrichment_results)) {
    all_enrichments[[contrast_name]] <- enrichment_results
    
    # Save individual results
    for (db_name in names(enrichment_results)) {
      write.csv(
        enrichment_results[[db_name]],
        file = paste0("pathway_analysis/", contrast_name, "_", db_name, ".csv"),
        row.names = FALSE
      )
    }
  }
}

# Combine all results
all_enrichments_df <- bind_rows(lapply(all_enrichments, function(x) {
  bind_rows(x)
}))

write.csv(all_enrichments_df, 
          "pathway_analysis/all_enrichments_combined.csv",
          row.names = FALSE)

cat("\n" , "=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("Over-representation analysis complete!\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")

# ==============================================================================
# Step 4: Visualize Enrichment Results
# ==============================================================================

cat("\nStep 4: Creating enrichment visualizations...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

# Dot plots for WikiPathways
pdf("figures/16_wikipathways_dotplots.pdf", width = 14, height = 20)

for (contrast_name in main_contrasts) {
  if (!is.null(all_enrichments[[contrast_name]][["WikiPathways"]])) {
    
    top_pathways <- all_enrichments[[contrast_name]][["WikiPathways"]] %>%
      arrange(p.adjust) %>%
      head(20) %>%
      mutate(
        GeneRatio_numeric = sapply(strsplit(GeneRatio, "/"), function(x) {
          as.numeric(x[1]) / as.numeric(x[2])
        }),
        Description_short = ifelse(nchar(Description) > 60,
                                   paste0(substr(Description, 1, 57), "..."),
                                   Description)
      )
    
    if (nrow(top_pathways) > 0) {
      p <- ggplot(top_pathways, 
                  aes(x = GeneRatio_numeric, y = reorder(Description_short, -p.adjust))) +
        geom_point(aes(size = Count, color = p.adjust)) +
        scale_color_gradient(low = "red", high = "blue") +
        theme_classic() +
        labs(
          title = paste("Top WikiPathways -", gsub("_", " ", contrast_name)),
          x = "Gene Ratio",
          y = "Pathway",
          color = "Adjusted\np-value",
          size = "Gene\nCount"
        ) +
        theme(
          axis.text.y = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, face = "bold")
        )
      
      print(p)
    }
  }
}

dev.off()

cat("WikiPathways dotplots saved\n")

# Dot plots for Hallmark
pdf("figures/17_hallmark_dotplots.pdf", width = 12, height = 16)

for (contrast_name in main_contrasts) {
  if (!is.null(all_enrichments[[contrast_name]][["Hallmark"]])) {
    
    hallmark_data <- all_enrichments[[contrast_name]][["Hallmark"]] %>%
      arrange(p.adjust) %>%
      head(20) %>%
      mutate(
        GeneRatio_numeric = sapply(strsplit(GeneRatio, "/"), function(x) {
          as.numeric(x[1]) / as.numeric(x[2])
        }),
        Description_clean = gsub("HALLMARK_", "", Description)
      )
    
    if (nrow(hallmark_data) > 0) {
      p <- ggplot(hallmark_data, 
                  aes(x = GeneRatio_numeric, y = reorder(Description_clean, -p.adjust))) +
        geom_point(aes(size = Count, color = p.adjust)) +
        scale_color_gradient(low = "red", high = "blue") +
        theme_classic() +
        labs(
          title = paste("Hallmark Gene Sets -", gsub("_", " ", contrast_name)),
          x = "Gene Ratio",
          y = "Gene Set",
          color = "Adjusted\np-value",
          size = "Gene\nCount"
        ) +
        theme(
          axis.text.y = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, face = "bold")
        )
      
      print(p)
    }
  }
}

dev.off()

cat("Hallmark dotplots saved\n")

# ==============================================================================
# Step 5: Compare Pathways Across Drugs
# ==============================================================================

cat("\nStep 5: Comparing pathways across drugs...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

drug_contrasts_noDHT <- c("Bica_vs_Control_noDHT", 
                          "Enza_vs_Control_noDHT", 
                          "Apalu_vs_Control_noDHT")

# Extract WikiPathways for comparison
wikipath_comparison <- lapply(drug_contrasts_noDHT, function(contrast) {
  if (!is.null(all_enrichments[[contrast]][["WikiPathways"]])) {
    all_enrichments[[contrast]][["WikiPathways"]] %>%
      select(Description, p.adjust, Count) %>%
      mutate(Drug = gsub("_vs_Control_noDHT", "", contrast))
  }
}) %>%
  bind_rows()

if (nrow(wikipath_comparison) > 0) {
  
  # Create pathway matrix
  pathway_matrix <- wikipath_comparison %>%
    select(Description, Drug, p.adjust) %>%
    pivot_wider(names_from = Drug, values_from = p.adjust, values_fill = 1) %>%
    column_to_rownames("Description")
  
  # Transform to -log10(FDR) for better visualization
  pathway_matrix_log <- -log10(pathway_matrix)
  
  # Keep only pathways present in at least 2 drugs
  pathway_sums <- rowSums(pathway_matrix < 0.05)
  pathway_matrix_filtered <- pathway_matrix_log[pathway_sums >= 2, ]
  
  if (nrow(pathway_matrix_filtered) > 0) {
    
    pdf("figures/18_pathway_comparison_heatmap.pdf", width = 10, height = 12)
    
    pheatmap(
      pathway_matrix_filtered,
      color = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      clustering_method = "complete",
      fontsize_row = 8,
      fontsize_col = 10,
      main = "Pathway Enrichment Across Drugs\n(-log10 adjusted p-value)",
      border_color = "grey60",
      cellwidth = 40,
      cellheight = 10
    )
    
    dev.off()
    
    cat("Pathway comparison heatmap saved\n")
    
    # Save comparison matrix
    write.csv(pathway_matrix_filtered, 
              "pathway_analysis/pathway_comparison_across_drugs.csv")
  }
}

# ==============================================================================
# Step 6: Classify Pathways
# ==============================================================================

cat("\nStep 6: Classifying pathways (shared vs drug-specific)...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

if (exists("wikipath_comparison") && nrow(wikipath_comparison) > 0) {
  
  pathway_classification <- wikipath_comparison %>%
    filter(p.adjust < 0.05) %>%
    group_by(Description) %>%
    summarize(
      Drugs = paste(sort(unique(Drug)), collapse = ", "),
      N_Drugs = n_distinct(Drug),
      Min_FDR = min(p.adjust),
      Total_Genes = sum(Count),
      .groups = "drop"
    ) %>%
    mutate(
      Category = case_when(
        N_Drugs == 3 ~ "Common to all drugs",
        N_Drugs == 2 ~ "Shared by two drugs",
        N_Drugs == 1 ~ "Drug-specific"
      )
    ) %>%
    arrange(Category, Min_FDR)
  
  write.csv(pathway_classification,
            "pathway_analysis/pathway_classification.csv",
            row.names = FALSE)
  
  cat("\nPathway classification:\n")
  print(table(pathway_classification$Category))
  
  # Visualize classification
  pdf("figures/19_pathway_classification_summary.pdf", width = 10, height = 6)
  
  classification_summary <- pathway_classification %>%
    count(Category) %>%
    mutate(Category = factor(Category, 
                             levels = c("Common to all drugs", 
                                        "Shared by two drugs", 
                                        "Drug-specific")))
  
  ggplot(classification_summary, aes(x = Category, y = n, fill = Category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = n), vjust = -0.5, size = 5) +
    theme_minimal() +
    labs(
      title = "Classification of Enriched Pathways",
      x = "Category",
      y = "Number of Pathways"
    ) +
    scale_fill_manual(values = c(
      "Common to all drugs" = "#E41A1C",
      "Shared by two drugs" = "#377EB8",
      "Drug-specific" = "#4DAF4A"
    )) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  
  dev.off()
  
  cat("Pathway classification plot saved\n")
}

# ==============================================================================
# Step 7: Pathway Network Analysis
# ==============================================================================

cat("\nStep 7: Creating pathway similarity network...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

if (exists("wikipath_comparison") && nrow(wikipath_comparison) > 0) {
  
  # Get all significant pathways
  sig_pathways <- wikipath_comparison %>%
    filter(p.adjust < 0.05) %>%
    pull(Description) %>%
    unique()
  
  if (length(sig_pathways) > 1) {
    
    # Calculate Jaccard similarity between pathways
    pathway_genes <- lapply(sig_pathways, function(pathway) {
      gene_string <- wikipathways_human %>%
        filter(gs_name == pathway) %>%
        pull(entrez_gene) %>%
        unique()
      return(gene_string)
    })
    names(pathway_genes) <- sig_pathways
    
    # Compute Jaccard similarity matrix
    jaccard_similarity <- function(set1, set2) {
      intersection <- length(intersect(set1, set2))
      union <- length(union(set1, set2))
      return(intersection / union)
    }
    
    n_pathways <- length(sig_pathways)
    similarity_matrix <- matrix(0, nrow = n_pathways, ncol = n_pathways)
    rownames(similarity_matrix) <- sig_pathways
    colnames(similarity_matrix) <- sig_pathways
    
    for (i in 1:n_pathways) {
      for (j in i:n_pathways) {
        if (i == j) {
          similarity_matrix[i, j] <- 1
        } else {
          sim <- jaccard_similarity(pathway_genes[[i]], pathway_genes[[j]])
          similarity_matrix[i, j] <- sim
          similarity_matrix[j, i] <- sim
        }
      }
    }
    
    # Create edge list (only keep similarities > 0.2)
    edge_list <- which(similarity_matrix > 0.2 & similarity_matrix < 1, arr.ind = TRUE)
    
    if (nrow(edge_list) > 0) {
      edges_df <- data.frame(
        Source = rownames(similarity_matrix)[edge_list[, 1]],
        Target = rownames(similarity_matrix)[edge_list[, 2]],
        Jaccard = similarity_matrix[edge_list],
        stringsAsFactors = FALSE
      )
      
      # Add gene overlap information
      edges_df$Shared_Genes <- sapply(1:nrow(edges_df), function(i) {
        length(intersect(
          pathway_genes[[edges_df$Source[i]]],
          pathway_genes[[edges_df$Target[i]]]
        ))
      })
      
      write.csv(edges_df, 
                "pathway_analysis/pathway_network_edges.csv",
                row.names = FALSE)
      
      # Create node attributes
      nodes_df <- data.frame(
        Pathway = sig_pathways,
        stringsAsFactors = FALSE
      )
      
      nodes_df <- nodes_df %>%
        left_join(
          wikipath_comparison %>%
            group_by(Description) %>%
            summarize(
              Drugs = paste(sort(unique(Drug)), collapse = ", "),
              N_Drugs = n_distinct(Drug),
              Min_FDR = min(p.adjust),
              .groups = "drop"
            ),
          by = c("Pathway" = "Description")
        )
      
      write.csv(nodes_df,
                "pathway_analysis/pathway_network_nodes.csv",
                row.names = FALSE)
      
      cat("Pathway network files created for Cytoscape\n")
      cat("  Nodes:", nrow(nodes_df), "\n")
      cat("  Edges:", nrow(edges_df), "\n")
    }
  }
}

# ==============================================================================
# Step 8: Gene Set Enrichment Analysis (GSEA)
# ==============================================================================

cat("\nStep 8: Running Gene Set Enrichment Analysis (GSEA)...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

run_GSEA_analysis <- function(DE_data, contrast_name) {
  
  cat("\nGSEA for:", contrast_name, "\n")
  
  # Create ranked gene list
  ranked_genes <- DE_data %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      rank_metric = sign(log2FoldChange) * (-log10(pvalue))
    ) %>%
    arrange(desc(rank_metric))
  
  # Convert to Entrez IDs
  gene_conversion <- convert_to_entrez(ranked_genes$GeneName)
  
  ranked_genes_entrez <- ranked_genes %>%
    inner_join(gene_conversion, by = c("GeneName" = "Symbol")) %>%
    distinct(EntrezID, .keep_all = TRUE)
  
  gene_list <- ranked_genes_entrez$rank_metric
  names(gene_list) <- ranked_genes_entrez$EntrezID
  
  cat("  Ranked genes:", length(gene_list), "\n")
  
  gsea_results <- list()
  
  # WikiPathways GSEA
  tryCatch({
    gsea_wiki <- GSEA(
      geneList = gene_list,
      TERM2GENE = data.frame(
        term = rep(names(wikipathways_list), sapply(wikipathways_list, length)),
        gene = unlist(wikipathways_list)
      ),
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      eps = 0
    )
    
    if (!is.null(gsea_wiki) && nrow(gsea_wiki@result) > 0) {
      gsea_results[["WikiPathways"]] <- gsea_wiki@result %>%
        mutate(Database = "WikiPathways", Contrast = contrast_name)
      cat("  WikiPathways GSEA:", nrow(gsea_wiki@result), "pathways\n")
    }
  }, error = function(e) {
    cat("  WikiPathways GSEA: Error -", e$message, "\n")
  })
  
  # KEGG GSEA
  tryCatch({
    gsea_kegg <- GSEA(
      geneList = gene_list,
      TERM2GENE = data.frame(
        term = rep(names(kegg_list), sapply(kegg_list, length)),
        gene = unlist(kegg_list)
      ),
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      eps = 0
    )
    
    if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
      gsea_results[["KEGG"]] <- gsea_kegg@result %>%
        mutate(Database = "KEGG", Contrast = contrast_name)
      cat("  KEGG GSEA:", nrow(gsea_kegg@result), "pathways\n")
    }
  }, error = function(e) {
    cat("  KEGG GSEA: Error -", e$message, "\n")
  })
  
  # Hallmark GSEA
  tryCatch({
    gsea_hallmark <- GSEA(
      geneList = gene_list,
      TERM2GENE = data.frame(
        term = rep(names(hallmark_list), sapply(hallmark_list, length)),
        gene = unlist(hallmark_list)
      ),
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      eps = 0
    )
    
    if (!is.null(gsea_hallmark) && nrow(gsea_hallmark@result) > 0) {
      gsea_results[["Hallmark"]] <- gsea_hallmark@result %>%
        mutate(Database = "Hallmark", Contrast = contrast_name)
      cat("  Hallmark GSEA:", nrow(gsea_hallmark@result), "gene sets\n")
    }
  }, error = function(e) {
    cat("  Hallmark GSEA: Error -", e$message, "\n")
  })
  
  return(gsea_results)
}

# Run GSEA for main contrasts
all_gsea_results <- list()

for (contrast_name in main_contrasts) {
  gsea_result <- run_GSEA_analysis(DE_results[[contrast_name]], contrast_name)
  
  if (length(gsea_result) > 0) {
    all_gsea_results[[contrast_name]] <- gsea_result
    
    # Save individual results
    for (db_name in names(gsea_result)) {
      write.csv(
        gsea_result[[db_name]],
        file = paste0("pathway_analysis/GSEA_", contrast_name, "_", db_name, ".csv"),
        row.names = FALSE
      )
    }
  }
}

# Combine all GSEA results
all_gsea_df <- bind_rows(lapply(all_gsea_results, function(x) {
  bind_rows(x)
}))

if (nrow(all_gsea_df) > 0) {
  write.csv(all_gsea_df,
            "pathway_analysis/all_GSEA_results.csv",
            row.names = FALSE)
}

cat("\n" , "=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("GSEA analysis complete!\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")

# ==============================================================================
# Step 9: Visualize GSEA Results
# ==============================================================================

cat("\nStep 9: Creating GSEA visualizations...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

if (nrow(all_gsea_df) > 0) {
  
  pdf("figures/20_GSEA_dotplots.pdf", width = 14, height = 18)
  
  for (contrast_name in main_contrasts) {
    if (!is.null(all_gsea_results[[contrast_name]][["Hallmark"]])) {
      
      hallmark_gsea <- all_gsea_results[[contrast_name]][["Hallmark"]] %>%
        arrange(p.adjust) %>%
        head(20) %>%
        mutate(
          Description_clean = gsub("HALLMARK_", "", Description),
          Description_clean = gsub("_", " ", Description_clean)
        )
      
      if (nrow(hallmark_gsea) > 0) {
        p <- ggplot(hallmark_gsea,
                    aes(x = NES, y = reorder(Description_clean, NES))) +
          geom_point(aes(size = setSize, color = p.adjust)) +
          geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
          scale_color_gradient(low = "red", high = "blue") +
          theme_classic() +
          labs(
            title = paste("GSEA: Hallmark Gene Sets -", gsub("_", " ", contrast_name)),
            x = "Normalized Enrichment Score (NES)",
            y = "Gene Set",
            color = "Adjusted\np-value",
            size = "Set Size"
          ) +
          theme(
            axis.text.y = element_text(size = 12),
            plot.title = element_text(hjust = 0.5, face = "bold")
          )
        
        print(p)
      }
    }
  }
  
  dev.off()
  
  cat("GSEA dotplots saved\n")
}

# ==============================================================================
# Step 10: Biological Interpretation Summary
# ==============================================================================

cat("\nStep 10: Generating biological interpretation...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

interpretation_text <- paste0(
  "=" %>% rep(78) %>% paste(collapse = ""),
  "\n",
  "BIOLOGICAL INTERPRETATION OF PATHWAY ANALYSIS\n",
  "=" %>% rep(78) %>% paste(collapse = ""),
  "\n\n"
)

if (exists("pathway_classification")) {
  
  # Common pathways
  common_pathways <- pathway_classification %>%
    filter(Category == "Common to all drugs") %>%
    arrange(Min_FDR) %>%
    head(10)
  
  if (nrow(common_pathways) > 0) {
    interpretation_text <- paste0(
      interpretation_text,
      "COMMON MECHANISMS (All Three Drugs):\n",
      "-" %>% rep(78) %>% paste(collapse = ""),
      "\n\n",
      "These pathways represent the core anti-androgen mechanism shared across\n",
      "all three drugs. They likely reflect fundamental responses to AR inhibition.\n\n",
      "Top Common Pathways:\n"
    )
    
    for (i in 1:min(10, nrow(common_pathways))) {
      interpretation_text <- paste0(
        interpretation_text,
        "  ", i, ". ", common_pathways$Description[i],
        " (FDR: ", formatC(common_pathways$Min_FDR[i], format = "e", digits = 2), ")\n"
      )
    }
    interpretation_text <- paste0(interpretation_text, "\n")
  }
  
  # Drug-specific pathways
  for (drug in c("Bica", "Enza", "Apalu")) {
    drug_specific <- pathway_classification %>%
      filter(Category == "Drug-specific", grepl(drug, Drugs)) %>%
      arrange(Min_FDR) %>%
      head(5)
    
    if (nrow(drug_specific) > 0) {
      drug_full <- case_when(
        drug == "Bica" ~ "Bicalutamide",
        drug == "Enza" ~ "Enzalutamide",
        drug == "Apalu" ~ "Apalutamide"
      )
      
      interpretation_text <- paste0(
        interpretation_text,
        drug_full, "-SPECIFIC MECHANISMS:\n",
        "-" %>% rep(78) %>% paste(collapse = ""),
        "\n\n"
      )
      
      for (i in 1:nrow(drug_specific)) {
        interpretation_text <- paste0(
          interpretation_text,
          "  ", i, ". ", drug_specific$Description[i],
          " (FDR: ", formatC(drug_specific$Min_FDR[i], format = "e", digits = 2), ")\n"
        )
      }
      interpretation_text <- paste0(interpretation_text, "\n")
    }
  }
}

interpretation_text <- paste0(
  interpretation_text,
  "=" %>% rep(78) %>% paste(collapse = ""),
  "\n",
  "RECOMMENDATIONS FOR FOLLOW-UP:\n",
  "=" %>% rep(78) %>% paste(collapse = ""),
  "\n\n",
  "1. Validate top common pathway genes with qRT-PCR or Western blot\n",
  "2. Investigate drug-specific mechanisms for combination therapy potential\n",
  "3. Examine genes in shared pathways for biomarker development\n",
  "4. Consider functional assays for pathways unique to most effective drug\n",
  "5. Explore network interactions using Cytoscape visualization\n\n"
)

cat(interpretation_text)
writeLines(interpretation_text, "pathway_analysis/BIOLOGICAL_INTERPRETATION.txt")

# ==============================================================================
# Step 11: Prepare Data for Cytoscape
# ==============================================================================

cat("\nStep 11: Preparing files for Cytoscape visualization...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

# Create gene-pathway associations for network
if (exists("all_enrichments_df") && nrow(all_enrichments_df) > 0) {
  
  gene_pathway_edges <- all_enrichments_df %>%
    filter(p.adjust < 0.05) %>%
    dplyr::select(Description, geneID, Contrast, Database, p.adjust) %>%
    separate_rows(geneID, sep = "/") %>%
    rename(Gene = geneID, Pathway = Description)
  
  write.csv(gene_pathway_edges,
            "pathway_analysis/cytoscape_gene_pathway_edges.csv",
            row.names = FALSE)
  
  pathway_nodes <- all_enrichments_df %>%
    filter(p.adjust < 0.05) %>%
    dplyr::select(Description, Database, Contrast, p.adjust, Count) %>%
    group_by(Description, Database) %>%
    summarize(
      Min_FDR = min(p.adjust),
      Max_Count = max(Count),
      N_Contrasts = n(),
      Contrasts = paste(unique(Contrast), collapse = ", "),
      .groups = "drop"
    )
  
  write.csv(pathway_nodes,
            "pathway_analysis/cytoscape_pathway_nodes.csv",
            row.names = FALSE)
  
  cat("Cytoscape files created\n")
  cat("  Gene-pathway edges:", nrow(gene_pathway_edges), "\n")
  cat("  Pathway nodes:", nrow(pathway_nodes), "\n")
}

# ==============================================================================
# Step 12: Save All Results
# ==============================================================================

cat("\nStep 12: Saving all pathway enrichment results...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

# Check which objects exist before saving
objects_to_save <- c("all_enrichments", "all_enrichments_df", 
                     "all_gsea_results", "all_gsea_df",
                     "DE_results", "DE_summary_combined")

# Add optional objects if they exist
if (exists("pathway_classification")) {
  objects_to_save <- c(objects_to_save, "pathway_classification")
  cat("  Including pathway_classification\n")
} else {
  cat("  Note: pathway_classification not created (no pathways to classify)\n")
}

if (exists("wikipath_comparison")) {
  objects_to_save <- c(objects_to_save, "wikipath_comparison")
  cat("  Including wikipath_comparison\n")
} else {
  cat("  Note: wikipath_comparison not created\n")
}

cat("Saving objects:", paste(objects_to_save, collapse = ", "), "\n")

# Save
save(list = objects_to_save,
     file = "module2_pathway_enrichment_results.RData")

cat("All results saved to: module2_pathway_enrichment_results.RData\n")

# ==============================================================================
# Final Summary
# ==============================================================================

summary_text <- paste0(
  "\n",
  "=" %>% rep(78) %>% paste(collapse = ""),
  "\n",
  "MODULE 2 COMPLETE - PATHWAY ENRICHMENT ANALYSIS\n",
  "=" %>% rep(78) %>% paste(collapse = ""),
  "\n\n",
  "METHODS:\n",
  "  - Over-Representation Analysis (ORA) with hypergeometric test\n",
  "  - Gene Set Enrichment Analysis (GSEA)\n",
  "  - Pathway classification (shared vs drug-specific)\n",
  "  - Network analysis (Jaccard similarity)\n\n",
  "DATABASES:\n",
  "  - WikiPathways\n",
  "  - KEGG\n",
  "  - Reactome\n",
  "  - Hallmark Gene Sets\n",
  "  - Gene Ontology (BP, MF, CC)\n\n",
  "KEY OUTPUTS:\n",
  "-" %>% rep(78) %>% paste(collapse = ""),
  "\n",
  "Results:\n",
  "  - all_enrichments_combined.csv\n",
  "  - all_GSEA_results.csv\n",
  #"  - pathway_classification.csv\n",
  "  - pathway_comparison_across_drugs.csv\n",
  "  - BIOLOGICAL_INTERPRETATION.txt\n",
  "  - Cytoscape network files\n\n",
  "Figures:\n",
  "  - 16_wikipathways_dotplots.pdf\n",
  "  - 17_hallmark_dotplots.pdf\n",
  "  - 18_pathway_comparison_heatmap.pdf\n",
  "  - 19_pathway_classification_summary.pdf\n",
  "  - 20_GSEA_dotplots.pdf\n\n",
  "NEXT STEP: Run Module 3 for network visualization (requires Cytoscape)\n",
  "=" %>% rep(78) %>% paste(collapse = ""),
  "\n"
)

cat(summary_text)
writeLines(summary_text, "pathway_analysis/MODULE2_SUMMARY.txt")

cat("\nModule 2 complete!\n\n")