# thyroiditis.R
# Copyright (C) 2026 Dougan Lab
# Author: Enoch Chung
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

# 1 Import 

## 1.1 import functions
source("EC_functions.R")
source("LA_Common_Plotting_Utilities_Edits.R")

## 1.2 Package Import
source("load_packages.R")

## 1.3 import data 
# Set Directory 
setwd('thyroid_data/')
dir <- getwd()

# Load counts matrix
counts <- read.csv("thyroid_counts.csv", row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
colnames(counts) <- gsub("\\.", "-", colnames(counts))

# Load metadata
metadata <- read.csv("thyroid_metadata.csv", row.names = 1,  check.names = FALSE)

# Ensure both contain the same cells
shared_cells <- intersect(colnames(counts), rownames(metadata))
counts <- counts[, shared_cells]
metadata <- metadata[shared_cells, ]

# Create Seurat object
lechner2023_f <- CreateSeuratObject(counts = counts, meta.data = metadata)
original_obj_lechner2023_f <- lechner2023_f

# 2 processing 
## 2.1 harmonization
lechner2023_f <- autoTune4(lechner2023_f)
lechner2023_f <- autoPilot(lechner2023_f, resolution = 0.5)

get_all_markers(lechner2023_f)

cd8_lechner <- subsetter(lechner2023_f, c(0, 1, 2, 3, 4, 5, 11, 12, 13, 14), original_obj_lechner2023_f)

cd8_lechner <- autoTune4(cd8_lechner)
cd8_lechner <- autoPilot(cd8_lechner, resolution = 0.5)

get_all_markers(cd8_lechner)

cd4_lechner <- subsetter(cd8_lechner, c(0, 1, 3, 4, 5, 7, 11, 12), original_obj_lechner2023_f)

cd4_lechner <- autoTune4(cd4_lechner)
cd4_lechner <- autoPilot(cd4_lechner, resolution = 0.5)

get_all_markers(cd4_lechner)

cd4.1_lechner <- subsetter(cd4_lechner, c(0, 1, 2, 3, 5, 6, 7, 8, 10), original_obj_lechner2023_f)

cd4.1_lechner <- autoTune4(cd4.1_lechner)
cd4.1_lechner <- autoPilot(cd4.1_lechner, resolution = 0.5)

get_all_markers(cd4.1_lechner)

cd4.2_lechner <- subsetter(cd4.1_lechner, c(0, 1, 2, 3, 4, 5, 7, 8, 9), original_obj_lechner2023_f)

cd4.2_lechner <- autoTune4(cd4.2_lechner)
cd4.2_lechner <- autoPilot(cd4.2_lechner, resolution = 0.8)

get_all_markers(cd4.2_lechner)

cd4.3_lechner <- cd8_subsetter(cd4.2_lechner, original_obj_lechner2023_f, threshold = 1.5)

cd4.3_lechner <- autoTune4(cd4.3_lechner)
cd4.3_lechner <- autoPilot(cd4.3_lechner, resolution = 0.4)

get_all_markers(cd4.3_lechner)

# 3 plotting 
seurat_obj <- cd4.3_lechner
cl <- as.character(seurat_obj$seurat_clusters)
cl[cl == "5"] <- "4"
seurat_obj$seurat_clusters <- factor(cl)
Idents(seurat_obj) <- seurat_obj$seurat_clusters
table(seurat_obj$seurat_clusters)

cluster_names <- c(
  "0" = "Tfh",
  "1" = "Naive CD4", 
  "2" = "Th1",
  "3" = "Resident memory",
  "4" = "Tph",
  "5" = ".",
  "6" = "Proflierating",
  "7" = "Treg",
  "8" = "CD8"
)
seurat_obj@meta.data$cluster_names <- cluster_names[as.character(seurat_obj$seurat_clusters)]

seurat_obj <- subset(seurat_obj, idents = "8", invert = TRUE)

## 3.1 IL-21 graph
cluster_colors_cd4 <- c("0" = "#ADACCB", "1" = "#F9989A" , "2" =  "#68397A", 
                        "4" = "#F5CB1A", "3" = "#CB449F", "5" = "mediumseagreen",
                        "6" = "palevioletred", "7" = "skyblue4", "8" = "brown",
                        "9" =  "skyblue2", "10" = "burlywood2", "11" = "#3D5EA9" , 
                        "12" = "khaki", "13" = "darkviolet", "14" = "turquoise", 
                        "15" = "#CBC8BB", "16" =  "royalblue"  , "17" = "#DF0000", 
                        "18" = "peru", "19" = "orange", "20" = "grey30",
                        "21" = "springgreen3" , "22" = "seagreen", "23" = "maroon4",
                        "24" = "black", "25" = "hotpink") 

il21_expr <- FetchData(seurat_obj, vars = c("IL21", "cluster_names"))

pct_df <- il21_expr %>%
  group_by(cluster_names) %>%
  summarise(
    pct_positive = mean(IL21 > 0) * 100
  ) %>%
  as.data.frame()

pct_df$cluster_names <- factor(pct_df$cluster_names, levels = cluster_names[as.character(0:17)])

ggplot(pct_df, aes(x = cluster_names, y = pct_positive, fill = cluster_names)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = setNames(cluster_colors_cd4[as.character(0:17)], cluster_names[as.character(0:17)])) +
  ggtitle("Percentage of IL-21+ Cells by Cluster") +
  xlab("Cluster") +
  ylab("% IL-21+ Cells") +
  scale_y_continuous(limits = c(0, 40), expand = c(0, 0)) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

## 3.2 heatmap
key_markers <- c("PDCD1", "CXCR5", "CXCL13", "IL21", "FOXP3", "TBX21", "MKI67", "GATA3", "IFNG", "CXCR3", "CCR6", "CCR7", "CD69")
avg_exp_scaled <- t(scale(t(AverageExpression(seurat_obj, features = key_markers, group.by = "seurat_clusters")$RNA)))

cluster_label_map <- setNames(
  seurat_obj@meta.data[["cluster_names"]],
  as.character(seurat_obj@meta.data[["seurat_clusters"]])
)
cluster_label_map <- cluster_label_map[!duplicated(names(cluster_label_map))]

cluster_label_map <- cluster_label_map[sub("^g", "", colnames(avg_exp_scaled))]

Heatmap(avg_exp_scaled,
        name = "Scaled\nExpression",
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        cluster_rows = FALSE, cluster_columns = FALSE,
        column_labels = cluster_label_map,
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 9))


## 3.3 UMAP
plot_umap_labeled <- function(seurat_obj, label_col = "cluster_names", title = "CD4+ T Cell Clusters") {
  
  if (!label_col %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Column '%s' not found in meta.data. Available columns: %s",
                 label_col, paste(colnames(seurat_obj@meta.data), collapse = ", ")))
  }
  
  umap_df <- as.data.frame(Embeddings(seurat_obj, "umap"))
  colnames(umap_df) <- c("UMAP_1", "UMAP_2")
  umap_df$label <- seurat_obj@meta.data[[label_col]]
  umap_df$cluster_id <- as.character(seurat_obj@meta.data[["seurat_clusters"]])
  
  label_to_color <- setNames(
    cluster_colors_cd4[umap_df$cluster_id],
    umap_df$label
  )
  label_to_color <- label_to_color[!duplicated(names(label_to_color))]
  
  centroids <- umap_df %>%
    group_by(label) %>%
    summarise(across(starts_with("UMAP"), median))
  
  ggplot(umap_df, aes(UMAP_1, UMAP_2, color = label)) +
    geom_point(size = 0.7, alpha = 0.9, stroke = 0) +
    geom_label(
      data = centroids,
      aes(x = UMAP_1, y = UMAP_2, label = label),
      color = "black",
      fill = "white",
      size = 3.2,
      linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = label_to_color) +
    labs(title = title, x = "UMAP_1", y = "UMAP_2", color = "Cell Type") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))
}

plot_umap_labeled(seurat_obj)

## 3.4 Expansion 
seurat_obj$condi <- case_when(
  seurat_obj$condition %in% "IRAE" ~ "IrAE",
  seurat_obj$condition %in% "Normal" ~ "Control",
  seurat_obj$condition %in% "HT" ~ "other",
  TRUE ~ NA_character_
)

pt_col        <- "sample"
condition_col <- "condi"
cluster_col   <- "cluster_names"
conditions_to_compare <- c("Control", "IrAE")
clusters_to_plot <- NULL                      
cond_fill  <- c("#F8C5C8", "#E63946")
cond_point <- c("#E89AA0", "#C81E2C")
names(cond_fill)  <- conditions_to_compare
names(cond_point) <- conditions_to_compare

md <- seurat_obj@meta.data %>%
  dplyr::select(all_of(c(pt_col, condition_col, cluster_col))) %>%
  dplyr::rename(pt = !!pt_col, condition = !!condition_col, cluster = !!cluster_col) %>%
  filter(condition %in% conditions_to_compare) %>%
  mutate(condition = factor(condition, levels = conditions_to_compare))

totals <- md %>%
  group_by(pt, condition) %>%
  summarise(total_cells = n(), .groups = "drop")

counts <- md %>%
  group_by(pt, condition, cluster) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  complete(cluster, nesting(pt, condition), fill = list(n_cells = 0))

plot_df <- counts %>%
  left_join(totals, by = c("pt", "condition")) %>%
  mutate(pct = 100 * n_cells / total_cells,
         condition = factor(condition, levels = conditions_to_compare))

ctrl_cond <- conditions_to_compare[1]
box_cond  <- conditions_to_compare[2]

plot_cluster <- function(cl) {
  df       <- filter(plot_df, cluster == cl)
  df_ctrl  <- filter(df, condition == ctrl_cond)
  df_box   <- filter(df, condition == box_cond)
  
  ggplot() +
    geom_crossbar(data = df_ctrl,
                  aes(x = condition, y = pct, ymin = pct, ymax = pct),
                  width = 0.6, color = "grey20", fatten = 1.2) +
    geom_jitter(data = df_ctrl,
                aes(x = condition, y = pct, color = condition),
                width = 0.1, size = 2.5, alpha = 0.85) +
    geom_boxplot(data = df_box,
                 aes(x = condition, y = pct, fill = condition),
                 outlier.shape = NA, alpha = 0.9, width = 0.6) +
    geom_jitter(data = df_box,
                aes(x = condition, y = pct, color = condition),
                width = 0.15, size = 2, alpha = 0.85) +
    scale_fill_manual(values = cond_fill) +
    scale_color_manual(values = cond_point) +
    scale_x_discrete(limits = conditions_to_compare, drop = FALSE) +
    scale_y_continuous(limits = c(0, 50), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Condition", y = "% of total cells (per sample)",
         title = paste("Cluster", cl)) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none")
}

cls <- if (is.null(clusters_to_plot)) sort(unique(plot_df$cluster)) else clusters_to_plot
plots <- lapply(cls, plot_cluster)
names(plots) <- cls

plots
