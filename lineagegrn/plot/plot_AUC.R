#!/usr/bin/env Rscript
# Load necessary libraries
library(furrr)
library(ggraph)
library(igraph)
library(ape)

# Define FindNodeDepth function
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

# Define fixed_color_scale
fixed_color_scale <- scale_color_gradientn(
    colors = c("#ffff99", "#ff7f10", "#822688"),
    limits = c(0.55, 0.95),
    oob = scales::squish
)

# Define plot_AUC_on_lineage function
plot_AUC_on_lineage <- function(tree_newick, auc_data, output_file, node_size = 6, edge_width = 0.3, total_time = 1.2, label_color = FALSE) {
    tree <- read.tree(text = tree_newick)
    tree_graph <- as.igraph(tree)

    nodes_time_df <- data.frame(node_id = 1:(tree$Nnode + length(tree$tip.label)), node_name = c(tree$tip.label, tree$node.label))
    times <- numeric()

    for (i in nodes_time_df$node_id) {
        time <- FindNodeDepth(tree, node = i, ncells = length(tree$tip.label))
        times <- c(times, time)
    }

    nodes_time <- times
    names(nodes_time) <- nodes_time_df$node_name

    V(tree_graph)$Time <- nodes_time[V(tree_graph)$name]
    V(tree_graph)$type <- V(tree_graph)$name
    V(tree_graph)$node_label <- V(tree_graph)$name

    colnames(auc_data) <- c("node_id", "value")
    auc_data$value <- scales::rescale(auc_data$value, to = c(min(auc_data$value), max(auc_data$value)))
    AUC <- auc_data$value
    names(AUC) <- auc_data$node_id

    graph_output <- ggraph(tree_graph, layout = "dendrogram", height = Time)
    graph_output <- graph_output + geom_edge_bend(aes(width = factor(0.4), fontface = "plain"))
    graph_output <- graph_output + geom_node_point(aes(color = AUC), size = node_size, shape = "circle")
    graph_output <- graph_output + geom_node_text(aes(label = node_label), nudge_x = 0, nudge_y = -0.1, size = 2)
    graph_output <- graph_output + geom_node_text(aes(label = as.character(round(auc_data$value, 3))), nudge_x = 0, nudge_y = 0, size = 2)
    graph_output <- graph_output + fixed_color_scale
    graph_output <- graph_output + ylim(c(total_time, 0)) + ylab("Time")
    graph_output <- graph_output + scale_edge_width_manual(values = edge_width, guide = "none")
    graph_output <- graph_output + scale_size_manual(values = 2, guide = "none")
    graph_output <- graph_output + theme(
        axis.line.y = element_line(linewidth = 0.1),
        axis.ticks.y = element_line(linewidth = 0.1),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        panel.background = element_rect(fill = NA, color = NA, linewidth = 0.1),
        axis.title.x = element_blank(),
        text = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6)
    )

    if (!label_color) {
        graph_output <- graph_output + theme_graph() + theme(legend.position = "none")
    }

    ggsave(output_file,
        plot = graph_output, width = 3.5, height = 2, dpi = 300,
        units = "in", limitsize = FALSE,
        bg = "transparent"
    )

    cat("Plot saved to:", output_file, "\n")
}

# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
    stop("Usage: Rscript plot_auc.R <tree_newick> <auc_file> <output_file> <label_color>", call. = FALSE)
}

# Parse arguments
tree_newick <- args[1]
auc_file <- args[2]
output_file <- args[3]
label_color <- as.logical(args[4])

# Read input data
auc_data <- read.csv(auc_file, row.names = NULL)

# Run the function
plot_AUC_on_lineage(tree_newick, auc_data, output_file, label_color = label_color)
