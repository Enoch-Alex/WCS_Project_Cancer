# ==============================================================================
# Module 3: Network Visualization with Cytoscape
# ==============================================================================
# Project: Comparing anti-androgen drug mechanisms in prostate cancer
# Goal: Create interactive network visualizations
# Requirements: Cytoscape must be installed and running
# ==============================================================================

library(RCy3)
library(igraph)
library(tidyverse)

# Load results from Module 2
load("module2_pathway_enrichment_results.RData")

cat("=" %>% rep(78) %>% paste(collapse = ""), "\n")
cat("MODULE 3: NETWORK VISUALIZATION WITH CYTOSCAPE\n")
cat("=" %>% rep(78) %>% paste(collapse = ""), "\n\n")

# ==============================================================================
# Step 1: Check Cytoscape Connection
# ==============================================================================

cat("\nStep 1: Checking Cytoscape connection...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

# Test connection
connection_ok <- FALSE
tryCatch({
  result <- cytoscapePing()
  cat("✓ Cytoscape is connected!\n")
  connection_ok <- TRUE
}, error = function(e) {
  cat("✗ Cannot connect to Cytoscape\n")
  cat("  Error:", e$message, "\n")
})

if (!connection_ok) {
  cat("\n")
  cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
  cat("ERROR: Cannot connect to Cytoscape\n")
  cat("=" %>% rep(60) %>% paste(collapse = ""), "\n\n")
  cat("Troubleshooting:\n")
  cat("  1. Make sure Cytoscape is running\n")
  cat("  2. Wait for 'CyREST API running' in bottom-left corner\n")
  cat("  3. Try restarting both Cytoscape and R\n\n")
  stop("Cytoscape connection failed")
}

# Try to get version info (but don't fail if it errors)
cat("\nCytoscape version information:\n")
tryCatch({
  version <- cytoscapeVersionInfo()
  
  # Handle both list and vector formats
  if (is.list(version)) {
    cat("  Cytoscape:", version$cytoscapeVersion, "\n")
    cat("  API:", version$apiVersion, "\n")
  } else {
    cat("  Version info:", version, "\n")
  }
}, error = function(e) {
  cat("  (Version info unavailable - this is OK)\n")
})

cat("\n✓ Ready to create networks!\n\n")

# ==============================================================================
# Step 2: Create Drug-Specific Gene Networks
# ==============================================================================

cat("\nStep 2: Creating gene networks for each drug...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

create_gene_network <- function(contrast_name, network_name) {
  
  cat("\nCreating network:", network_name, "\n")
  
  # Get top DE genes
  de_genes <- DE_results[[contrast_name]] %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    head(100)
  
  if (nrow(de_genes) < 10) {
    cat("  Not enough significant genes. Skipping.\n")
    return(NULL)
  }
  
  cat("  Using", nrow(de_genes), "genes\n")
  
  # Check if gene-pathway file exists
  if (!file.exists("pathway_analysis/cytoscape_gene_pathway_edges.csv")) {
    cat("  No pathway association file found. Skipping.\n")
    return(NULL)
  }
  
  # Load gene-pathway associations
  gene_pathway_raw <- read.csv("pathway_analysis/cytoscape_gene_pathway_edges.csv")
  
  cat("  Loaded", nrow(gene_pathway_raw), "gene-pathway associations\n")
  
  # Convert Entrez IDs to Gene Symbols
  cat("  Converting Entrez IDs to gene symbols...\n")
  
  entrez_to_symbol <- mapIds(org.Hs.eg.db,
                             keys = as.character(unique(gene_pathway_raw$Gene)),  # ✅ FIXED
                             column = "SYMBOL",
                             keytype = "ENTREZID",
                             multiVals = "first")
  
  # Create conversion dataframe
  conversion_df <- data.frame(
    EntrezID = names(entrez_to_symbol),
    GeneSymbol = as.character(entrez_to_symbol),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(GeneSymbol))
  
  cat("  Converted", nrow(conversion_df), "IDs to symbols\n")
  
  # Add gene symbols to pathway data
  gene_pathway <- gene_pathway_raw %>%
    mutate(EntrezID = as.character(Gene)) %>%
    left_join(conversion_df, by = "EntrezID") %>%
    filter(!is.na(GeneSymbol)) %>%
    rename(Gene_Original = Gene, Gene = GeneSymbol)
  
  # Filter for this contrast and genes in our DE results
  gene_pathway <- gene_pathway %>%
    filter(Contrast == contrast_name, Gene %in% de_genes$GeneName)
  
  if (nrow(gene_pathway) == 0) {
    cat("  No pathway associations found after ID conversion. Skipping.\n")
    return(NULL)
  }
  
  cat("  Found", nrow(gene_pathway), "associations for", 
      length(unique(gene_pathway$Gene)), "genes\n")
  
  # Create edges: genes that share pathways
  gene_pairs <- gene_pathway %>%
    dplyr::select(Gene, Pathway) %>%
    inner_join(gene_pathway %>% dplyr::select(Gene, Pathway), 
               by = "Pathway",
               relationship = "many-to-many") %>%
    filter(Gene.x < Gene.y) %>%
    group_by(Gene.x, Gene.y) %>%
    summarize(Shared_Pathways = n(), .groups = "drop") %>%
    filter(Shared_Pathways >= 2)
  
  if (nrow(gene_pairs) == 0) {
    cat("  No gene connections found (no genes share >=2 pathways). Skipping.\n")
    return(NULL)
  }
  
  cat("  Creating", nrow(gene_pairs), "edges\n")
  
  # Prepare node attributes
  nodes_df <- data.frame(
    id = de_genes$GeneName,
    GeneName = de_genes$GeneName,
    log2FC = de_genes$log2FoldChange,
    padj = de_genes$padj,
    negLog10padj = -log10(de_genes$padj),
    Regulation = de_genes$Regulation,
    stringsAsFactors = FALSE
  )
  
  # Only keep nodes that have edges
  nodes_in_network <- unique(c(gene_pairs$Gene.x, gene_pairs$Gene.y))
  nodes_df <- nodes_df %>%
    filter(id %in% nodes_in_network)
  
  cat("  Network has", nrow(nodes_df), "nodes\n")
  
  # Prepare edge attributes
  edges_df <- data.frame(
    source = gene_pairs$Gene.x,
    target = gene_pairs$Gene.y,
    interaction = "shares_pathway",
    Shared_Pathways = gene_pairs$Shared_Pathways,
    stringsAsFactors = FALSE
  )
  
  # Create network in Cytoscape
  tryCatch({
    net_suid <- createNetworkFromDataFrames(
      nodes = nodes_df,
      edges = edges_df,
      title = network_name,
      collection = "Drug Networks"
    )
    
    # Style the network
    # Node color by regulation
    setNodeColorMapping(
      "Regulation",
      mapping.type = "discrete",
      table.column.values = c("UP", "DOWN", "NS"),
      colors = c("#E41A1C", "#377EB8", "#CCCCCC"),
      network = net_suid
    )
    
    # Node size by significance
    setNodeSizeMapping(
      "negLog10padj",
      mapping.type = "continuous",
      table.column.values = c(min(nodes_df$negLog10padj), 
                              max(nodes_df$negLog10padj)),
      sizes = c(30, 80),
      network = net_suid
    )
    
    # Edge width by shared pathways
    if (max(edges_df$Shared_Pathways) > min(edges_df$Shared_Pathways)) {
      setEdgeLineWidthMapping(
        "Shared_Pathways",
        mapping.type = "continuous",
        table.column.values = c(min(edges_df$Shared_Pathways), 
                                max(edges_df$Shared_Pathways)),
        widths = c(1, 5),
        network = net_suid
      )
    }
    
    # Apply layout
    layoutNetwork("force-directed", network = net_suid)
    
    cat("  ✓ Network created successfully\n")
    
    # Export image
    output_file <- paste0("figures/21_network_", gsub(" ", "_", network_name), ".png")
    exportImage(output_file, type = "PNG", resolution = 300, network = net_suid)
    cat("  ✓ Image exported\n")
    
    return(network_name)
    
  }, error = function(e) {
    cat("  ✗ Error creating network:", e$message, "\n")
    return(NULL)
  })
}



# Make sure Cytoscape is running
cytoscapePing()

# Clear and redefine the function above, then:
drug_networks <- list()

drug_networks[["Bicalutamide_DHT"]] <- create_gene_network(
  "Bica_vs_Control_DHT",
  "Bicalutamide_DHT"
)

drug_networks[["Enzalutamide_DHT"]] <- create_gene_network(
  "Enza_vs_Control_DHT",
  "Enzalutamide_DHT"
)

drug_networks[["Apalutamide_DHT"]] <- create_gene_network(
  "Apalu_vs_Control_DHT",
  "Apalutamide_DHT"
)

# Save session
saveSession("results/Cytoscape_Drug_Networks.cys")

cat("\n✓ Networks complete!\n")

# ==============================================================================
# Step 3: Save Cytoscape Session
# ==============================================================================

cat("\nStep 3: Saving Cytoscape session...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

session_file <- "results/Cytoscape_Complete_Session.cys"

tryCatch({
  saveSession(session_file)
  cat("  ✓ Session saved to:", session_file, "\n")
  cat("  You can reopen this file in Cytoscape later\n")
}, error = function(e) {
  cat("  ✗ Could not save session:", e$message, "\n")
})

# ==============================================================================
# Step 4: Network Statistics
# ==============================================================================

cat("\nStep 4: Computing network statistics...\n")
cat("-" %>% rep(60) %>% paste(collapse = ""), "\n")

tryCatch({
  network_names <- getNetworkList()
  
  if (length(network_names) > 0) {
    cat("Created", length(network_names), "networks:\n")
    for (net in network_names) {
      cat("  -", net, "\n")
    }
  } else {
    cat("  No networks were created\n")
  }
}, error = function(e) {
  cat("  Could not list networks:", e$message, "\n")
})

# ==============================================================================
# Final Summary
# ==============================================================================

summary_text <- paste0(
  "\n",
  "=" %>% rep(78) %>% paste(collapse = ""),
  "\n",
  "MODULE 3 COMPLETE - NETWORK VISUALIZATION\n",
  "=" %>% rep(78) %>% paste(collapse = ""),
  "\n\n",
  "NETWORKS CREATED:\n",
  "-" %>% rep(78) %>% paste(collapse = ""),
  "\n",
  "  - Bicalutamide_DHT\n",
  "  - Enzalutamide_DHT\n",
  "  - Apalutamide_DHT\n\n",
  "FILES CREATED:\n",
  "-" %>% rep(78) %>% paste(collapse = ""),
  "\n",
  "Networks:\n",
  "  - figures/21_network_Bicalutamide_DHT.png\n",
  "  - figures/21_network_Enzalutamide_DHT.png\n",
  "  - figures/21_network_Apalutamide_DHT.png\n\n",
  "Session:\n",
  "  - results/Cytoscape_Complete_Session.cys\n\n",
  "NEXT STEPS:\n",
  "-" %>% rep(78) %>% paste(collapse = ""),
  "\n",
  "1. Open Cytoscape session file to explore networks interactively\n",
  "2. Use Cytoscape apps for additional analysis\n",
  "3. Export high-resolution images for publication\n",
  "4. Review pathway enrichment results from Module 2\n\n",
  "=" %>% rep(78) %>% paste(collapse = ""),
  "\n"
)

cat(summary_text)
writeLines(summary_text, "results/MODULE3_SUMMARY.txt")

cat("\n✓ All modules complete! 🎉\n\n")