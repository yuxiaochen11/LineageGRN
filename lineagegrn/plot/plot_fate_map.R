#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(ggraph)
  library(ggplot2)
  library(igraph)
})


FindNodeDepth <- function(tree, node, ncells) {
  root_node <- ncells + 1
  tree_edges <- tree$edge
  frame <- cbind(tree$edge, tree$edge.length)
  depth <- 0
  i <- node
  while (i != (ncells + 1)) {
    depth <- depth + frame[frame[, 2] == i, 3]
    i <- tree_edges[tree_edges[, 2] == i, ][1]
  }
  return(depth)
}

NODES_COLORS <- c(
  "#E78585", "#C8FFC2", "#F3E985", "#85C4F1", "#B1A867", "#D66B6B", "#FFD966", "#406E92", "#FFF2CC", "#B9F380", "#6F87D1", "#AF5454", "#8BD353", "#7982CB", "#3D6A85",
  "#247BA0", "#96C777", "#D1C488", "#FCE883", "#8EA4F5", "#EEDD82", "#C5B489", "#A4956A", "#6DB4D6", "#3D87A1", "#9E8C55", "#5F95B1", "#7AB6D6", "#3D6A85", "#6497B1",
  "#A4E78D", "#C25959", "#9F3E3E", "#9CB2F3", "#B4CBF9"
)

plot_fate_map <- function(tree_newick, developmental_time, node_size, edge_width, output_file) {
  tree <- read.tree(text = tree_newick)
  tree_graph <- as.igraph(tree)
  nodes_num <- 2 * length(tree$tip.label) + 1
  nodes_color <- NODES_COLORS[1:nodes_num]
  nodes_time_df <- data.frame(
    node_id = 1:(tree$Nnode + length(tree$tip.label)),
    node_name = c(tree$tip.label, tree$node.label)
  )
  times <- c()
  for (i in nodes_time_df$node_id) {
    time <- FindNodeDepth(tree, node = i, ncells = length(tree$tip.label))
    times <- c(times, time)
  }
  nodes_time <- times
  names(nodes_time) <- nodes_time_df$node_name

  V(tree_graph)$Time <- nodes_time[V(tree_graph)$name]
  V(tree_graph)$type <- V(tree_graph)$name
  V(tree_graph)$node_label <- V(tree_graph)$name
  V(tree_graph)$Cat <- "Tip"
  V(tree_graph)$Cat[V(tree_graph)$name %in% tree$node.label] <- "Node"

  graph_out <- ggraph(tree_graph, layout = "dendrogram", height = Time) +
    geom_edge_bend(
      aes(width = factor(1), fontface = "plain"),
      # arrow = arrow(type = "closed", length = unit(2, "pt"), angle = 45),
      start_cap = circle(5, "pt"),
      end_cap = circle(5, "pt")
    ) +
    geom_node_point(aes(color = type, shape = Cat), size = node_size) +
    geom_node_text(aes(label = node_label), nudge_x = 0, nudge_y = -0.08) +
    scale_shape_manual(values = c(18, 15)) +
    scale_color_manual(values = nodes_color) +
    ylim(c(developmental_time, 0)) +
    ylab("Time") +
    scale_edge_width_manual(values = edge_width, guide = "none") +
    theme(
      panel.background = element_rect(fill = NA, color = NA, size = 10),
      legend.position = "none",
      axis.title.x = element_blank(),
      text = element_text(size = 12)
    )


  ggsave(output_file, plot = graph_out, width = 12, height = 8, dpi = 300, units = "in")
  cat("Plot saved to:", output_file, "\n")
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  print(length(args))
  stop("Usage: Rscript plot_fate_map.R <tree_newick> <developmental_time> <node_size> <edge_width> <output_file>")
}

tree_newick <- args[1]
developmental_time <- as.numeric(args[2])
node_size <- as.numeric(args[3])
edge_width <- as.numeric(args[4])
output_file <- args[5]



plot_fate_map(tree_newick, developmental_time, node_size, edge_width, output_file)
