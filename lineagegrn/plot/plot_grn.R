#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(umap)
  library(ggplot2)
  library(ggraph)
  library(igraph)
  library(tidyverse)
  library(RColorBrewer)
  library(readr)
  library(ggrepel)
  library(factoextra)
})


# Function for UMAP visualization and clustering
gene_umap <- function(counts_file, target_gene_file, regulator_file) {
  set.seed(123)
  target_gene <- read.csv(target_gene_file)[[1]]
  regulator <- read.csv(regulator_file)[[1]]
  genes <- c(target_gene, regulator)
  expression_counts_matrix_ <- read_csv(counts_file, col_names = FALSE)
  expression_counts_matrix <- expression_counts_matrix_[, 3:ncol(expression_counts_matrix_)]
  expression_counts_matrix <- as.data.frame(lapply(expression_counts_matrix, function(x) as.numeric(as.character(x))))

  umap_config <- umap.defaults
  umap_config$n_neighbors <- 3
  umap_config$min_dist <- 0.5
  umap_result <- umap(expression_counts_matrix, config = umap_config)
  umap_df <- as.data.frame(umap_result$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$Gene <- genes

  k <- 4
  kmeans_result <- kmeans(umap_df[, c("UMAP1", "UMAP2")], centers = k)
  umap_df$Cluster <- as.factor(kmeans_result$cluster)

  ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster, label = Gene)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    ggtitle("UMAP Visualization of Genes with Clustering") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_text_repel(max.overlaps = 20, size = 2, box.padding = 0.5, point.padding = 0.5)

  return(umap_df)
}

# Function for plotting the gene network
network_plot <- function(gene_interactions_file, umap_df, output_pdf_file) {
  gene_interactions <- read_csv(gene_interactions_file)
  network_data <- data.frame(from = gene_interactions[, 2], to = gene_interactions[, 3])
  umap_data <- data.frame(Gene = umap_df$Gene, UMAP1 = umap_df$UMAP1, UMAP2 = umap_df$UMAP2)

  graph <- graph_from_data_frame(d = network_data, vertices = umap_data, directed = TRUE)

  degree_graph <- degree(graph)
  for (i in 1:length(degree_graph)) {
    if (startsWith(names(degree_graph)[i], "r")) {
      degree_graph[i] <- degree_graph[i] / 1.5
    }
  }

  V(graph)$UMAP1 <- umap_df$UMAP1
  V(graph)$UMAP2 <- umap_df$UMAP2
  V(graph)$Size <- 40
  V(graph)$Color <- degree_graph

  p <- ggraph(graph, layout = "manual", x = V(graph)$UMAP1, y = V(graph)$UMAP2) +
    geom_edge_fan(alpha = 0.01, color = "gray") +
    geom_node_point(aes(size = V(graph)$Size, color = Color), alpha = 1, shape = 16) +
    labs(color = "Degree") +
    guides(size = "none") +
    scale_size_continuous(range = c(1, 5)) +
    scale_color_viridis_c() +
    theme_void()

  pdf(output_pdf_file, width = 5, height = 4)
  print(p)
  dev.off()
}

# Main function for command-line arguments
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    stop("Usage: Rscript gene_network_analysis.R <counts_file> <target_gene_file> <regulator_file> <gene_interactions_file> <output_pdf_file>")
  }

  counts_file <- args[1]
  target_gene_file <- args[2]
  regulator_file <- args[3]
  gene_interactions_file <- args[4]
  output_pdf_file <- args[5]

  # Read genes from the provided file

  # Generate UMAP data
  umap_df <- gene_umap(counts_file, target_gene_file, regulator_file)

  # Generate network plot
  network_plot(gene_interactions_file, umap_df, output_pdf_file)
}

# Run the script
if (interactive() == FALSE) {
  main()
}
