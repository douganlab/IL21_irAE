# LA_Common_Plotting_Utilities.R
# Copyright (C) 2025 Lestat Ali
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

library(ggplot2)
library(ggrepel)

PREFERRED_PLOT_FORMAT <- "pdf"

EMBEDDING_WIDTH <- 6
EMBEDDING_HEIGHT <- 5.5

default_cluster_color_panel <- c("0" = "#ADACCB", "1" = "#F9989A" , "2" =  "#68397A", 
                                 "4" = "#F5CB1A", "3" = "#CB449F", "5" = "mediumseagreen",
                                 "6" = "palevioletred", "7" = "skyblue4", "8" = "brown",
                                 "9" =  "skyblue2", "10" = "burlywood2", "11" = "#3D5EA9" , 
                                 "12" = "khaki", "13" = "darkviolet", "14" = "turquoise", 
                                 "15" = "#CBC8BB", "16" =  "royalblue"  , "17" = "#DF0000", 
                                 "18" = "peru", "19" = "orange", "20" = "grey30",
                                 "21" = "springgreen3" , "22" = "seagreen", "23" = "maroon4",
                                 "24" = "black", "25" = "hotpink") 
default_cluster_color_panel <- c(default_cluster_color_panel, rep("black", 25))
names(default_cluster_color_panel)[27:51] <- as.character(26:50)

default_category_color_panel <- c("0" = "firebrick2", "1" = "lightskyblue2", "2" =  "khaki2", 
                                  "3" = "dodgerblue4", "4" = "plum3", "5" = "darkslategray",
                                  "6" = "peachpuff", "7" = "chocolate1", "8" = "lightpink", 
                                  "9" = "antiquewhite3", "10" = "steelblue3",
                                  "11" = "darkseagreen", "12" = "dimgray")
default_category_color_panel <- c(default_category_color_panel, paste0("grey", seq(5, 100, 10)), rep("black", 20))
default_category_color_panel <- unname(default_category_color_panel)

.save_plot <- function(filename, plot, width, height) {
  if (PREFERRED_PLOT_FORMAT == "pdf") {
    filename <- paste0(filename, ".pdf")
    ggsave(filename, plot, device = cairo_pdf, width = width, height = height)
  } else if (PREFERRED_PLOT_FORMAT == "png") {
    filename <- paste0(filename, ".png")
    ggsave(filename, plot, device = "png", width = width, height = height, dpi = 300, type = "cairo")
  } else {
    print(paste0("Error! Unknown file format:", PREFERRED_PLOT_FORMAT))
  }
}


theme_la <- function() {
  return(theme_bw() + theme(text = element_text(family = "Helvetica Neue"), 
                            legend.title = element_text(face = "plain"),
                            plot.title = element_text(face = "bold"), 
                            panel.grid = element_blank(),
                            panel.border = element_blank(),
                            axis.title.x = element_text(size = 10),
                            axis.title.y = element_text(size = 10),
                            axis.text.x = element_text(size = 7),
                            axis.text.y = element_text(size = 7),
                            strip.text = element_text(size = 12, colour = "white", face = "bold", family = "Helvetica Neue"),
                            strip.background = element_rect(fill = "grey15")))
}


theme_la_bar <- function() {
  return(theme_la() + theme(axis.line = element_line(size = 0.25),
                            axis.text.y = element_text(size = 11),
                            axis.text.x = element_text(size = 11),
                            axis.title.y = element_text(size = 13, margin = margin(r = 12))))
}


theme_la_umap <- function() {
  return(theme_la() + theme(axis.line = element_blank(),
                            axis.ticks = element_blank(),
                            axis.text.x = element_blank(),
                            axis.text.y = element_blank(),
                            axis.title.x = element_blank(),
                            axis.title.y = element_blank(),
                            panel.border = element_rect(fill = NA)))
}


draw_cluster_umap <- function(seurat_object, label = T, cluster_labels = NULL, label_sizes = 11, colors = NULL, pt_size = 0.3, alpha = 0.8, save = NULL) {
  .draw_cluster_2d(seurat_object, "umap", label, cluster_labels, label_sizes, colors, pt_size, alpha, save)
}

draw_cluster_tsne <- function(seurat_object, label = T, cluster_labels = NULL, label_sizes = 11, colors = NULL, pt_size = 0.3, alpha = 0.8, save = NULL) {
  .draw_cluster_2d(seurat_object, "tsne", label, cluster_labels, label_sizes, colors, pt_size, alpha, save)
}

.draw_cluster_2d <- function(seurat_object, reduction, label = T, cluster_labels = NULL, label_sizes = 11, colors = NULL, pt_size = 0.3, alpha = 0.8, save = NULL) {
  embeddings_df <- .get_embeddings(seurat_object, reduction)
  
  varx <- colnames(embeddings_df)[1]
  vary <- colnames(embeddings_df)[2]
  
  embeddings_df$Cluster <- Idents(seurat_object)
  embedding_centers_df <- aggregate(x = cbind(eval(parse(text = varx)), eval(parse(text = vary))) ~ Cluster, data = embeddings_df, FUN = mean)
  colnames(embedding_centers_df)[2:3] <- c(varx, vary)
  embedding_centers_df$LabelSize <- label_sizes
  if (is.null(cluster_labels)) {
    embedding_centers_df$ClusterLabel <- embedding_centers_df$Cluster
  } else {
    embedding_centers_df$ClusterLabel <- cluster_labels[embedding_centers_df$Cluster]
  }
  
  cluster_colors <- default_cluster_color_panel
  if (!is.null(colors)) {
    cluster_colors <- colors
  }
  
  plot <- ggplot(data = embeddings_df, aes(x = get(varx), y = get(vary))) +
    geom_point(shape = 16, alpha = alpha, size = pt_size, aes(color = Cluster, fill = Cluster)) +
    scale_color_manual(values = cluster_colors) + 
    scale_fill_manual(values = cluster_colors) +
    guides(color = "none", fill = "none") + 
    xlab(paste0(toupper(reduction), " Dim 1")) + ylab(paste0(toupper(reduction), " Dim 2")) +
    theme_la_umap()
  
  if (label) {
    embedding_centers_sub_df <- subset(embedding_centers_df, !is.na(LabelSize) & !is.na(ClusterLabel))
      plot <- plot + geom_point(data = embedding_centers_sub_df, mapping = aes(size = LabelSize), color = "grey15", alpha = 0.65) +
        geom_text(data = embedding_centers_sub_df, aes(label = ClusterLabel), color = "white", fontface = "bold") +
        #scale_size_area(max_size = max(label_sizes), guide = "none")
       scale_radius(range = c(5, 5 + max(label_sizes, na.rm = T)), guide = "none")
  }
  
  if (!is.null(save)) {
    filename <- paste0("plots/", seurat_object@project.name, "_", save, "_cluster_", reduction)
    .save_plot(filename, plot, width = EMBEDDING_WIDTH, height = EMBEDDING_HEIGHT)
  } else {
    plot
  }
}


draw_subset_cluster_umap <- function(seurat_object, cell_subset, label = T, colors = NULL, pt_size = 0.3, save = NULL) {
  .draw_subset_cluster_2d(seurat_object, "umap", cell_subset, label, colors, pt_size, save)
}

draw_subset_cluster_tsne <- function(seurat_object, cell_subset, label = T, colors = NULL, pt_size = 0.3, save = NULL) {
  .draw_subset_cluster_2d(seurat_object, "tsne", cell_subset, label, colors, pt_size, save)
}

.draw_subset_cluster_2d <- function(seurat_object, reduction, cell_subset, label = T, colors = NULL, pt_size = 0.3, save = NULL) {
  embeddings_df <- .get_embeddings(seurat_object, reduction)
  
  varx <- colnames(embeddings_df)[1]
  vary <- colnames(embeddings_df)[2]
  
  embeddings_df$Cluster <- Idents(seurat_object)
  embedding_centers_df <- aggregate(x = cbind(eval(parse(text = varx)), eval(parse(text = vary))) ~ Cluster, data = embeddings_df, FUN = mean)
  colnames(embedding_centers_df)[2:3] <- c(varx, vary)
  
  cluster_colors <- default_cluster_color_panel
  if (!is.null(colors)) {
    cluster_colors <- colors
  }
  
  shown_embeddings_df <- embeddings_df[rownames(embeddings_df) %in% cell_subset,]
  hidden_embeddings_df <- embeddings_df[!(rownames(embeddings_df) %in% cell_subset),]
  
  plot <- ggplot(data = embeddings_df, aes(x = get(varx), y = get(vary), color = Cluster, fill = Cluster)) +
    geom_point(data = hidden_embeddings_df, shape = 16, alpha = 0, color = NA, fill = NA, size = pt_size) +
    geom_point(data = shown_embeddings_df, shape = 16, alpha = 0.85, size = pt_size) +
    scale_color_manual(values = cluster_colors) + 
    scale_fill_manual(values = cluster_colors) +
    guides(color = "none", fill = "none") + 
    xlab(paste0(toupper(reduction), " Dim 1")) + ylab(paste0(toupper(reduction), " Dim 2")) +
    theme_la_umap()
  
  if (label) {
    plot <- plot + geom_text(data = embedding_centers_df, aes(label = Cluster), color = "black", fontface = "bold")
  }
  
  if (!is.null(save)) {
    filename <- paste0("plots/", seurat_object@project.name, "_", save, "_cluster_", reduction)
    .save_plot(filename, plot, width = EMBEDDING_WIDTH, height = EMBEDDING_HEIGHT)
  } else {
    plot
  }
}


draw_scaled_counts_umap <- function(seurat_object, genes, weights = NULL, assay = NULL, limit_pos = NULL, limit_neg = NULL, pt_size = 0.45, show_title = T, show_legend = T, save = NULL) {
  .draw_scaled_counts_2d(seurat_object, "umap", genes, weights, assay, limit_pos, limit_neg, pt_size, show_title, show_legend, save)
}

draw_scaled_counts_tsne <- function(seurat_object, genes, weights = NULL, assay = NULL, limit_pos = NULL, limit_neg = NULL, pt_size = 0.3, show_title = T, show_legend = T, save = NULL) {
  .draw_scaled_counts_2d(seurat_object, "tsne", genes, weights, assay, limit_pos, limit_neg, pt_size, show_title, show_legend, save)
}

.draw_scaled_counts_2d <- function(seurat_object, reduction, genes, weights = NULL, assay = NULL, limit_pos = NULL, limit_neg = NULL, pt_size = 0.3, show_title = T, show_legend = T, save = NULL) {
  reduction <- "umap"
  
  embeddings_df <- .get_embeddings(seurat_object, reduction)
  embeddings_df$Cluster <- Idents(seurat_object)
  embeddings_df$Score <- score_cells_for_markers(seurat_object, genes, weights, assay)
  
  min_score <- min(embeddings_df$Score)
  max_score <- max(embeddings_df$Score)
  
  if (!is.null(limit_pos)) {
    max_score <- limit_pos #min(limit_pos, max_score)
  }   
  
  if (!is.null(limit_neg)) {
    min_score <- limit_neg #max(limit_neg, min_score)
  } 
  
  mid_color <- ifelse(abs(min_score) < 0.1, "grey95", "white")
  
  concat_genes <- paste(genes, collapse = ", ")
  if (length(genes) > 3) {
    concat_genes <- paste0(paste(genes[1:3], collapse = ", "), "...")
  }
  
  title <- ifelse(is.null(save), concat_genes, save)
  if (!show_title) {
    title <- NULL
  }
  
  varx <- colnames(embeddings_df)[1]
  vary <- colnames(embeddings_df)[2]
  
  print(summary(embeddings_df$Score))
  
  plot <- ggplot(mapping = aes(x = get(varx), y = get(vary))) +
    geom_point(data = subset(embeddings_df, abs(Score) < 0.5), shape = 16, alpha = 0.5, size = 0.8 * pt_size, color = "grey95") +
    geom_point(data = subset(embeddings_df, !(abs(Score) < 0.5)), mapping = aes(color = Score), shape = 16, alpha = 1, size = pt_size) +
    scale_color_gradient2(name = "Expression", low = "navy", high = "red3", mid = mid_color, guide = "none",
                          midpoint = 0, limits = c(0.9 * min_score, 0.9 * max_score), oob = scales::squish) +
    ggtitle(title) +
    xlab(paste0(toupper(reduction), " Dim 1")) + ylab(paste0(toupper(reduction), " Dim 2")) +
    theme_la_umap() +
    theme(title = element_text(size = 11, face = "bold"))
          # panel.border = element_blank(),
          # axis.title.x = element_blank(),
          # axis.title.y = element_blank(),
          # legend.title = element_text(size = 8),
          # legend.text = element_text(size = 7))
  
  if (show_legend)  {
    plot <- plot + guides(color = guide_colorbar(ticks.colour = "black", frame.colour = "black"))
  }
  
  if (!is.null(save)) {
    filename <- paste0("plots/", seurat_object@project.name, "_", save, "_scaled", "_", reduction)
    .save_plot(filename, plot, width = EMBEDDING_WIDTH, height = EMBEDDING_HEIGHT)
  } else {
    plot
  }
}


score_cells_for_markers <- function(seurat_object, markers, weights = NULL, assay = NULL) {
  invalid_markers <- setdiff(markers, rownames(seurat_object))
  if (length(invalid_markers) > 0) {
    print(paste0("Ignoring ", length(invalid_markers), " unmeasured features."))
  }
  
  valid_markers <- intersect(markers, rownames(seurat_object))
  if (length(valid_markers) == 0) {
    stop("No valid markers for plotting!")
  } 
  
  if (!is.null(weights) & length(valid_markers) != length(weights)) {
    stop("Markers and weights need to be of the same length.")
  }
  
  assay_used <- ifelse(is.null(assay), "RNA", assay)
  slot_used <- ifelse(grepl("SCT", assay_used), "scale.data", "data")
  is_scaled <- grepl("SCT", assay_used)
  
  marker_exp_mat <- as.matrix(GetAssayData(seurat_object, layer = slot_used, assay = assay_used)[valid_markers,])
  if (length(valid_markers) > 1) {
    marker_exp_mat <- t(marker_exp_mat)
  }
  
  zero_genes <- colSums(marker_exp_mat) == 0
  if (length(which(zero_genes)) > 0) {
    print(paste0("Ignoring zero-expression genes: ", colnames(marker_exp_mat)[zero_genes]))
    marker_exp_mat <- marker_exp_mat[, !zero_genes]
  }
  
  if (is_scaled) {
    scaled_marker_exp_mat <- marker_exp_mat
  } else {
    scaled_marker_exp_mat <- scale(marker_exp_mat)
  }
  
  if (!is.null(weights)) {
    valid_weights <- weights[which(markers %in% valid_markers)]
    scaled_marker_exp_mat <- t(t(scaled_marker_exp_mat) * valid_weights)
  }
  
  scores <- rowSums(scaled_marker_exp_mat)
  return(scores)
}


draw_highlight_umap <- function(seurat_object, highlight_cells, pt_size = 0.3, fill = "grey15", color = "grey65", alpha = 0.9, save = NULL) {
  .draw_highlight_2d(seurat_object, "umap", highlight_cells, pt_size, fill, color, alpha, save)
}

draw_highlight_tsne <- function(seurat_object, highlight_cells, pt_size = 0.3, fill = "grey15", color = "grey65", alpha = 0.9, save = NULL) {
  .draw_highlight_2d(seurat_object, "tsne", highlight_cells, pt_size, fill, color, alpha, save)
}

.draw_highlight_2d <- function(seurat_object, reduction, highlight_cells, pt_size = 0.3, fill = "grey15", color = "grey65", alpha = 0.9, save = NULL) {
  embeddings_df <- .get_embeddings(seurat_object, reduction)
  fg_embeddings_df <- embeddings_df[rownames(embeddings_df) %in% highlight_cells,]
  bg_embeddings_df <- embeddings_df[!(rownames(embeddings_df) %in% highlight_cells),]
  
  varx <- colnames(embeddings_df)[1]
  vary <- colnames(embeddings_df)[2]
  
  plot <- ggplot(data = embeddings_df, aes(x = get(varx), y = get(vary))) +
    geom_point(data = bg_embeddings_df, shape = 21, fill = "grey85", color = "grey85", alpha = 0.3, size = pt_size) +
    geom_point(data = fg_embeddings_df, shape = 21, fill = fill, color = color, alpha = alpha, size = pt_size * 1.25) +
    guides(colour = "none") +
    xlab(paste0(toupper(reduction), " Dim 1")) + ylab(paste0(toupper(reduction), " Dim 2")) +
    theme_la_umap()
  
  if (!is.null(save)) {
    filename <- paste0("plots/", seurat_object@project.name, "_", save, "_highlight_", reduction)
    .save_plot(filename, plot, width = EMBEDDING_WIDTH, height = EMBEDDING_HEIGHT)
  } else {
    plot
  }
}


draw_meta_cat_umap <- function(seurat_object, cat_meta_var, colors = NULL, pt_size = 0.3, alpha = 0.75, show_legend = T, legend_title = "", bg_fg = F, fg_multiplier = 1, save = NULL) {
  .draw_meta_cat_2d(seurat_object, "umap", cat_meta_var, colors, pt_size, alpha, show_legend, legend_title, bg_fg, fg_multiplier, save)
}

draw_meta_cat_tsne <- function(seurat_object, cat_meta_var, colors = NULL, pt_size = 0.3, alpha = 0.75, show_legend = T, legend_title = "", bg_df = F, fg_multiplier = 1, save = NULL) {
  .draw_meta_cat_2d(seurat_object, "tsne", cat_meta_var, colors, pt_size, alpha, show_legend, legend_title, bg_fg, fg_multiplier, save)
}

.draw_meta_cat_2d <- function(seurat_object, reduction, cat_meta_var, colors = NULL, pt_size = 0.3, alpha = 0.75, show_legend = T, legend_title = "", bg_fg = F, fg_multiplier = 1, save = NULL) {
  embeddings_df <- .get_embeddings(seurat_object, reduction)
  embeddings_df$Cluster <- Idents(seurat_object)
  embeddings_df$Category <- as.character(seurat_object@meta.data[, cat_meta_var])
  
  cat_colors <- default_category_color_panel
  if (!is.null(colors)) {
    cat_colors <- colors
  }
  
  varx <- colnames(embeddings_df)[1]
  vary <- colnames(embeddings_df)[2]
  
  if (bg_fg) {
   bg_categories <- c("None", "Unmatched", "Ambiguous", "Unknown", "No", "FALSE", "LO")
  } else {
    bg_categories <- c()
  }
  
  plot <- ggplot(data = embeddings_df, aes(x = get(varx), y = get(vary))) +
    geom_point(data = subset(embeddings_df, Category %in% bg_categories), shape = 16, aes(color = Category, fill = Category), alpha = alpha, size = pt_size) +
    geom_point(data = subset(embeddings_df, !(Category %in% bg_categories)), shape = 16, aes(color = Category, fill = Category), alpha = fg_multiplier * alpha, size = fg_multiplier * pt_size) +
    scale_color_manual(values = cat_colors, name = legend_title) + 
    scale_fill_manual(values = cat_colors, name = legend_title) +
    xlab(paste0(toupper(reduction), " Dim 1")) + ylab(paste0(toupper(reduction), " Dim 2")) +
    theme_la_umap() + 
    theme(legend.position = "bottom")
  
  if (show_legend) {
    plot <- plot + guides(color = guide_legend(override.aes = list(size = 3)), fill = guide_legend())
  } else {
    plot <- plot + guides(color = "none", fill = "none")
  }
  
  if (!is.null(save)) {
    filename <- paste0("plots/", seurat_object@project.name, "_", save, "_", cat_meta_var, "_", reduction)
    height <- ifelse(show_legend, EMBEDDING_HEIGHT + 1, EMBEDDING_HEIGHT)
    .save_plot(filename, plot, width = EMBEDDING_WIDTH, height = height)
  } else {
    plot
  }
}

draw_clonal_expansion_umap <- function(seurat_object, clonality_df, pt_size = 0.3, max_cells = 250, show_legend = T, save = NULL) {
  .draw_clonal_expansion_2d(seurat_object, "umap", clonality_df, pt_size, max_cells, show_legend, save)
}

draw_clonal_expansion_tsne <- function(seurat_object, clonality_df, pt_size = 0.3, max_cells = 250, show_legend = T, save = NULL) {
  .draw_clonal_expansion_2d(seurat_object, "tsne", clonality_df, pt_size, max_cells, show_legend, save)
}

.draw_clonal_expansion_2d <- function(seurat_object, reduction, clonality_df, pt_size = 0.3, max_cells = 250, show_legend = T, save = NULL) {
  embeddings_df <- .get_embeddings(seurat_object, reduction)
  embeddings_df$CloneFrequency <- clonality_df[rownames(embeddings_df), "CloneFrequency"]
  
  varx <- colnames(embeddings_df)[1]
  vary <- colnames(embeddings_df)[2]
  
  max_score <- max_cells
  
  plot <- ggplot(data = embeddings_df, aes(x = get(varx), y = get(vary))) +
    geom_point(data = subset(embeddings_df, is.na(CloneFrequency) | CloneFrequency == 0), color = "grey90", shape = 16, alpha = 0.95, size = pt_size) +
    geom_point(data = subset(embeddings_df, CloneFrequency == 1), color = "darkseagreen2", shape = 16, alpha = 0.95, size = pt_size) +
    geom_point(data = subset(embeddings_df, CloneFrequency > 1), aes(color = CloneFrequency), shape = 16, alpha = 0.95, size = pt_size) +
    scale_color_gradient(low = "lightsteelblue2", high = "steelblue4", name = "Clones",
                         na.value = "black", oob = scales::squish, limits = c(2, max_score)) +
    xlab(paste0(toupper(reduction), " Dim 1")) + ylab(paste0(toupper(reduction), " Dim 2")) +
    theme_la_umap() +
    theme(legend.position = c(0.9, 0.22))
  
  if (show_legend) {
    plot <- plot + guides(color = guide_colorbar(ticks.colour = "black", frame.colour = "black"))
  } else {
    plot <- plot + guides(color = guide_none())
  }
  
  if (!is.null(save)) {
    filename <- paste0("plots/", seurat_object@project.name, "_", save, "_clonality_", reduction)
    .save_plot(filename, plot, width = EMBEDDING_WIDTH, height = EMBEDDING_HEIGHT)
  } else {
    plot
  }
}

  
draw_volcano <- function(deg_df, maxp = NULL, minp = NULL, maxfc = NULL, title = "", 
                         label = T, labeled_up_genes = NULL, labeled_dn_genes = NULL, n_label = 10, label_size = 3,
                         up_color = "red2", dn_color = "dodgerblue3", save = NULL) {
  if (!is.data.frame(deg_df)) {
    stop("Error, first argument must be a data frame.")
  }
  
  most_sig_p <- min(deg_df$p_val_adj[deg_df$p_val_adj != 0], na.rm = T)
  deg_df$p_val_adj[deg_df$p_val_adj == 0] <- 0.9 * most_sig_p
  deg_df$Gene <- rownames(deg_df)
  sig_deg_df <- subset(deg_df, p_val_adj <= 0.05)
  up_sig_deg_df <- subset(sig_deg_df, avg_log2FC > 0)
  dn_sig_deg_df <- subset(sig_deg_df, avg_log2FC < 0)
  if (!is.null(labeled_up_genes)) {
    most_up_genes <- labeled_up_genes
  } else {
    most_up_genes <- select_important_genes(up_sig_deg_df, n_label) 
  }
  if (!is.null(labeled_dn_genes)) {
    most_dn_genes <- labeled_dn_genes
  } else {
    most_dn_genes <- select_important_genes(dn_sig_deg_df, n_label) 
  }
  up_non_deg_df <- subset(deg_df, p_val_adj > 0.05 & avg_log2FC >= 0)
  dn_non_deg_df <- subset(deg_df, p_val_adj > 0.05 & avg_log2FC < 0)
  
  
  if (is.null(maxp)) {
    maxp <- 1.02 * max(-log10(deg_df$p_val_adj), na.rm = T)
  }
  
  if (is.null(minp)) {
    minp <- - 0.01 * maxp
  }
  
  if (is.null(maxfc)) {
    maxfc <- 1.02 * max(abs(deg_df$avg_log2FC), na.rm = T)
  }
  
  volcano_plot <- ggplot(deg_df, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_jitter(data = up_non_deg_df, shape = 19, color = up_color, alpha = 0.3, size = 0.7, height = 0.1, width = 0.01) +
    geom_jitter(data = dn_non_deg_df, shape = 19, color = dn_color, alpha = 0.3, size = 0.7, height = 0.1, width = 0.01) +
    geom_point(data = up_sig_deg_df, shape = 19, color = up_color, alpha = 0.8, size = 1.2) +
    geom_point(data = dn_sig_deg_df, shape = 19, color = dn_color, alpha = 0.8, size = 1.2) +
    annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(0.05), fill = "snow", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    xlab("Fold Change (log2)") + ylab("P (-log10)") +
    ggtitle(title) +
    scale_x_continuous(expand = expansion(0), limits = c(-1 * maxfc, maxfc)) +
    scale_y_continuous(limits = c(minp, maxp)) +
    theme_bw() +
    theme(text = element_text(size = 12), #14
          axis.title = element_text(size = 12), #14
          plot.title = element_text(family = "Avenir Next", size = 20, face = "bold", hjust = 0),
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
    
  if (label) {
    volcano_plot <- volcano_plot + 
      geom_text_repel(data = subset(deg_df, Gene %in% most_up_genes), max.overlaps = 20, force = 20,
                      max.time = 2, force_pull = 0.5, segment.alpha = 0.3, fontface = "plain",
                       aes(label = Gene),  min.segment.length = 0.3, nudge_x = 0 + 0.7, nudge_y = 0.1, size = label_size) + #size 3.2

      geom_text_repel(data = subset(deg_df, Gene %in% most_dn_genes), max.overlaps = 10, force = 20,
                      max.time = 2, force_pull = 0.5, segment.alpha = 0.3, fontface = "plain",
                       aes(label = Gene),  nudge_y = -.1, nudge_x = -0.7, min.segment.length = 0.5, size = label_size)
  }
    
  if (!is.null(save)) {
    .save_plot(paste0("plots/", save, "_deg_volcano"), volcano_plot, width = 5.1, height = 5)
  } else {
    volcano_plot
  }
}


select_important_genes <- function(deg_df, n = 10) {
  if (nrow(deg_df) < 1) {
    return(NA)
  }
  
  # take only top 10% of genes by significance
  p_val_cutoff <- quantile(deg_df$p_val_adj, 0.25)
  deg_df <- subset(deg_df, p_val_adj <= p_val_cutoff)
  
  if (deg_df$avg_log2FC[1] > 0) {
    ordered_lfc_genes <- rownames(deg_df)[order(deg_df$avg_log2FC, decreasing = T)]
  } else {
    ordered_lfc_genes <- rownames(deg_df)[order(deg_df$avg_log2FC, decreasing = F)]
  }
  
  if (!grepl(ordered_lfc_genes[1], "_")) {
    ordered_lfc_genes <- ordered_lfc_genes[!startsWith(ordered_lfc_genes, "RP") &
                                             !startsWith(ordered_lfc_genes, "IG") &
                                             !startsWith(ordered_lfc_genes, "TRA") &
                                             !startsWith(ordered_lfc_genes, "TRB") &
                                             !startsWith(ordered_lfc_genes, "MT") &
                                             !startsWith(ordered_lfc_genes, "Rp") &
                                             !startsWith(ordered_lfc_genes, "mt")]
  }
  
  top_lfc_genes <- ordered_lfc_genes[1:(n+2)]
  
  ordered_pval_genes <- rownames(deg_df)[order(deg_df$p_val_adj)]
  
  if (!grepl(ordered_pval_genes[1], "_")) {
    ordered_pval_genes <- ordered_pval_genes[!startsWith(ordered_pval_genes, "RP") &
                                               !startsWith(ordered_pval_genes, "IG") &
                                               !startsWith(ordered_pval_genes, "TRA") &
                                               !startsWith(ordered_pval_genes, "TRB") &
                                               !startsWith(ordered_pval_genes, "MT") &
                                               !startsWith(ordered_pval_genes, "Rp") &
                                               !startsWith(ordered_pval_genes, "mt")]
  }
  
  top_pval_genes <- ordered_pval_genes[1:n]
  
  return(union(top_lfc_genes, top_pval_genes))
}


draw_violin_sandwich_rev <- function(seurat_object, genes, ordered_clusters, cluster_name_map, labels = NULL, colors = NULL, save = NULL) {
  count_data <- GetAssayData(seurat_object$RNA)
  cluster_ids <- Idents(seurat_object)
  gene_expr <- count_data[genes,]
  plottable_expr <- reshape2::melt(t(as.matrix(gene_expr)))
  plottable_expr$Cluster <- rep(cluster_ids, length(genes))
  colnames(plottable_expr) <- c("CellID", "Gene", "Expression", "ClusterID")
  plottable_expr$ClusterName <- cluster_name_map[plottable_expr$ClusterID]
  plottable_expr <- subset(plottable_expr, ClusterID %in% ordered_clusters)
  plottable_expr$ClusterID <- factor(plottable_expr$ClusterID, levels = ordered_clusters)
  
  
  if (is.null(labels)) {
    labels <- ordered_clusters
  }
  
  get.labels <- function(x) {
    return(as.data.frame(labels))
  }
  
  if (is.null(colors)) {
    colors <- default_cluster_color_panel
  }
  
  expr_violinp <- ggplot(data = plottable_expr, aes(x = Gene, y = Expression)) +
    annotate(geom = "segment", x = 0, xend = 0, y = 0, yend = 4) +
    geom_violin(aes(fill = ClusterID, color = ClusterID), size = 0.8, alpha = 0.9, width = 0.8, scale = "width", trim = T) + 
    facet_wrap(facets = vars(ClusterID), labeller = get.labels, ncol = 1, dir = "v", strip.position = "right") + 
    scale_colour_manual(values = colors) +
    scale_fill_manual(values = colors) +
    guides(fill = "none", color = "none") +
    scale_x_discrete(position = "top", expand = expansion(add = c(0, 0.5))) +
    scale_y_continuous(breaks = c(0, 2, 4)) +
    theme_la() + 
    theme(panel.spacing = unit(0, "cm"),
          plot.title = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank(),
          axis.text.x = element_text(size = 12, color = "black", face = "bold", family = "Helvetica Neue"),
          axis.text.y = element_text(size = 12),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  

  for (i in 1:length(genes)) {
    expr_violinp <- expr_violinp + geom_vline(xintercept = i + 0.5, color = "grey70", linetype = "dotted")
  }

  if (!is.null(save)) {
    .save_plot(paste0("plots/violins_", save), 
               plot = expr_violinp, 
               width = 1.1 * length(genes) / 1.55 ,
               height = 1.1 * length(ordered_clusters) / 1.55)
  } else {
    expr_violinp
  }
}


draw_gsea_pathway_running_plot <- function(pathway, deg_df, col = "black", save = NULL)  {
  
  rankings <- deg_df$avg_log2FC
  names(rankings) <- rownames(deg_df)
  
  rnk <- rank(-rankings)
  statsAdj <- rankings[order(rnk)]
  statsAdj <- sign(statsAdj) * (abs(statsAdj))
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- sort(unname(as.vector(na.omit(match(pathway, names(statsAdj))))))
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, returnAllExtremes = T)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  max <- max(c(abs(bottoms), abs(tops)), na.rm = T)
  max <- max(0.5, max)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, length(statsAdj) + 1), y = c(0, ys, 0))
  segments_df <- data.frame("Rank" = pathway, "Stat" = statsAdj[pathway])
  lfc_rank_df <- data.frame("Position" = rnk, "LFC" = rankings[names(rnk)])
  midpoint <- length(which(rankings > 0))
  diff <- (max(tops) - min(bottoms))/6
  x = y = NULL
  plot <- ggplot(toPlot, aes(x = x, y = y)) +
    geom_line(color = col, size = 0.75) +
    geom_hline(yintercept = 0, colour = "black") + 
    geom_segment(data = segments_df, 
                 mapping = aes(x = Rank, y = -diff, xend = Rank, yend = diff, alpha = abs(Stat)),
                 size = 0.5, colour = "black") +
    guides(alpha = "none") +
    scale_y_continuous(limits = c(-max, max)) +
    scale_x_continuous(expand = expansion(add = 25)) +
    geom_segment(data = lfc_rank_df, y = -Inf, yend = -max, 
                 aes(x = Position, xend = Position, colour = LFC), size = 1) +
    scale_colour_gradient2(limits = c(-0.5, 0.5),
                           oob = scales::squish, low = "dodgerblue3", high = "red2", mid = "white",
                           midpoint = 0, guide = "none") +
    theme_la() +
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = NA, color = "black"))
  if (!is.null(save)) {
    .save_plot(paste0("plots/", save, "_gsea"), plot, width = 3.2, height = 2)
  }
  
  return(plot)
}


draw_enrichment_cols <- function(gsea_df, save = NULL) {
  pathway_fontfaces <- rep("plain", nrow(gsea_df))
  pathway_fontfaces[gsea_df$Sig] <- "bold"
  pathway_fontcolors <- rep("black", nrow(gsea_df))
  pathway_fontcolors[!gsea_df$Sig] <- "grey80"
  gsea_col_plot <- ggplot(gsea_df, aes(y = pathway, x = NES)) +
    geom_col(aes(fill = Sig)) +
    geom_vline(xintercept = 0) +
    scale_fill_manual(values = c("grey90", "grey25"), guide = "none") +
    theme_la_bar() +
    theme(axis.text.y = element_text(face = pathway_fontfaces, color = pathway_fontcolors, size = 7),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank())
  if (!is.null(save)) {
    .save_plot(paste0("plots/", save, "_gsea_plot"), gsea_col_plot, width = 6, height = 6)
  } else {
    gsea_col_plot
  }
}


.get_embeddings <- function(seurat_object, reduction_prefix) {
  desired_reduction <- NULL
  if (reduction_prefix == "umap") {
    valid_reductions <- Reductions(seurat_object)
    desired_reduction <- "umap.rpca"
    if (!(desired_reduction %in% valid_reductions)) {
      desired_reduction <- "umap"
    }
  } else if (reduction_prefix == "tnse") {
    desired_reduction <- "tsne"
  } 
  
  if (is.null(desired_reduction)) {
    stop("Error! Invalid reduction prefix!")
  } else {
    return(as.data.frame(Embeddings(seurat_object, desired_reduction)))
  }
}

sem <- function(x, na.rm = TRUE) {
  if (na.rm) {
    return(sd(x, na.rm = T) / sqrt(length(na.omit(x))))
  } else {
    return(sd(x, na.rm = F) / sqrt(length(x)))
  }
}

