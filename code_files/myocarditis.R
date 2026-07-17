# myocarditis.R
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

## 1.3 data import 
import_h5_files <- function(dir_path, tissue_filter = "heart") {
  h5_files <- list.files(dir_path, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)
  
  if (length(h5_files) == 0) {
    message("No .h5 files found in: ", dir_path)
    return(NULL)
  }
  
  # Parse metadata from filenames
  metadata <- data.frame(
    file_path  = h5_files,
    gsm_id     = gsub("^(GSM\\d+)_.*", "\\1", basename(h5_files)),
    patient_id = gsub("^GSM\\d+_SIC_(\\d+)_.*", "\\1", basename(h5_files)),
    sample_id  = gsub("^GSM\\d+_SIC_\\d+_(\\d+)_.*", "\\1", basename(h5_files)),
    tissue     = gsub("^.*_(heart|lung|blood)_.*", "\\1", basename(h5_files))
  )
  
  # Filter to tissue of interest
  metadata <- metadata[metadata$tissue == tissue_filter, ]
  
  if (nrow(metadata) == 0) {
    message("No files found for tissue: ", tissue_filter)
    return(NULL)
  }
  
  message("Found ", nrow(metadata), " ", tissue_filter, " samples")
  
  seurat_list <- lapply(seq_len(nrow(metadata)), function(i) {
    message("Importing: ", basename(metadata$file_path[i]))
    counts <- Read10X_h5(metadata$file_path[i])
    seu <- CreateSeuratObject(counts = counts, project = metadata$gsm_id[i])
    seu$gsm_id     <- metadata$gsm_id[i]
    seu$patient <- metadata$patient_id[i]
    seu$sample  <- metadata$sample_id[i]
    seu$tissue     <- metadata$tissue[i]
    return(seu)
  })
  
  names(seurat_list) <- metadata$gsm_id
  return(seurat_list)
}
blum2024_list <- import_h5_files("GSE228597_RAW", tissue_filter = "heart")

merge_union(blum2024_list, "blum2024_f")

# 2 Processing
blum2024_f <- autoTune4(blum2024_f)
blum2024_f <- autoPilot(blum2024_f, resolution = 0.5)
get_all_markers(blum2024_f)

T_blum2024 <- subsetter(blum2024_f, c(5, 7, 9, 18), original_obj_blum2024_f)

T_blum2024 <- autoTune4(T_blum2024)
T_blum2024 <- autoPilot(T_blum2024, resolution = 0.5)
get_all_markers(T_blum2024)

cd4_blum2024 <- subsetter(T_blum2024, c(1, 3, 5, 8, 14), original_obj_blum2024_f)

cd4_blum2024 <- autoTune4(cd4_blum2024)
cd4_blum2024 <- autoPilot(cd4_blum2024, resolution = 0.5)
get_all_markers(cd4_blum2024)

cd4.1_blum2024 <- cd8_subsetter(cd4_blum2024, original_obj_blum2024_f, threshold = 2.5)

cd4.1_blum2024 <- autoTune4(cd4.1_blum2024)
cd4.1_blum2024 <- autoPilot(cd4.1_blum2024, resolution = 0.25)
get_all_markers(cd4.1_blum2024)

# 3 Plotting
seurat_obj <- cd4.1_blum2024

cluster_names <- c(
  "0" = "Resting CD4",
  "1" = "Proliferating", 
  "2" = "Memory CD4",
  "3" = "Activated",
  "4" = "Mixed"
)
seurat_obj@meta.data$cluster_names <- cluster_names[as.character(seurat_obj$seurat_clusters)]

seurat_obj <- subset(seurat_obj, idents = "4", invert = TRUE)

# IL-21 graph 
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

# Heatmap
key_markers <- c("PDCD1", "CXCR5", "CXCL13", "IL21", "IL23R", "IL17A", "RORC", "FOXP3", "TBX21", "MKI67", "GATA3", "IFNG", "CXCR3", "CCR6", "CCR7", "CD69")
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

# UMAP
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
































