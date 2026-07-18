# colitis_thomas.R
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

## 1.3 Data import

# Helper function to extract sample ID from folder name
extract_sample_id_from_folder <- function(folder_name) {
  # Look for SIC_XX or MC_X pattern
  if (grepl("SIC_\\d+|MC_\\d+", folder_name)) {
    match <- regmatches(folder_name, regexpr("SIC_\\d+|MC_\\d+", folder_name))
    return(match)
  }
  return(NA)
}

# Main function to load raw data
load_raw_data <- function(data_path = "data/colon_data2024") {
  
  raw_data_path <- file.path(data_path, "GSE206299_RAW")
  
  if (!dir.exists(raw_data_path)) {
    stop(paste("RAW data directory not found:", raw_data_path))
  }
  
  # Find all channel directories
  channel_dirs <- list.dirs(raw_data_path, recursive = FALSE)
  channel_dirs <- channel_dirs[grepl("channels", basename(channel_dirs))]
  
  cat("Found", length(channel_dirs), "channel directories\n")
  
  # Initialize list to store Seurat objects
  seurat_list <- list()
  
  for (channel_dir in channel_dirs) {
    
    # Get the sample subdirectories
    sample_subdirs <- list.dirs(channel_dir, recursive = FALSE)
    
    for (sample_subdir in sample_subdirs) {
      
      sample_name <- basename(sample_subdir)
      cat("Processing:", sample_name, "\n")
      
      # Extract sample ID from folder name
      sample_id <- extract_sample_id_from_folder(sample_name)
      
      # Skip if not a GEX sample
      if (is.na(sample_id) || !grepl("GEX", sample_name)) {
        cat("  Skipping (not GEX)\n")
        next
      }
      
      # Check for filtered or raw 10X files (prefer filtered)
      filtered_h5 <- file.path(sample_subdir, "filtered_feature_bc_matrix.h5")
      raw_h5 <- file.path(sample_subdir, "raw_feature_bc_matrix.h5")
      
      if (file.exists(filtered_h5)) {
        h5_file <- filtered_h5
        cat("  Using filtered matrix\n")
      } else if (file.exists(raw_h5)) {
        h5_file <- raw_h5
        cat("  Using raw matrix\n")
      } else {
        cat("  Skipping (no H5 file found)\n")
        next
      }
      
      tryCatch({
        # Load 10X data
        cat("  Loading 10X data...\n")
        counts <- Read10X_h5(h5_file)
        
        # Create Seurat object
        seurat_obj <- CreateSeuratObject(counts = counts, 
                                         project = sample_id, 
                                         min.cells = 3,
                                         min.features = 200)
        
        # Add sample ID
        seurat_obj$sample <- sample_id
        
        cat("  Loaded:", ncol(seurat_obj), "cells,", nrow(seurat_obj), "genes\n")
        
        # Create unique key for this sample
        # If multiple folders exist for same sample, append suffix
        key <- sample_id
        suffix <- 1
        while (key %in% names(seurat_list)) {
          suffix <- suffix + 1
          key <- paste0(sample_id, "_", suffix)
        }
        
        seurat_list[[key]] <- seurat_obj
        
      }, error = function(e) {
        cat("  Error loading", sample_name, ":", conditionMessage(e), "\n")
      })
    }
  }
  
  if (length(seurat_list) == 0) {
    stop("No valid samples were loaded!")
  }
  
  cat("\nLoaded", length(seurat_list), "sample(s):\n")
  for (name in names(seurat_list)) {
    cat("  ", name, ":", ncol(seurat_list[[name]]), "cells\n")
  }
  
  return(seurat_list)
}
thomas_2024_list <- load_raw_data("colon_data2024")

## 1.4 merging
merge_union(thomas_2024_list, "thomas2024_f")

# 2 Processing 
thomas2024_f <- autoTune4(thomas2024_f)
thomas2024_f <- autoPilot(thomas2024_f, resolution = 0.5)
get_all_markers(thomas2024_f)

cd8_thomas <- subsetter(thomas2024_f, c(0, 1, 3, 4, 8, 21), original_obj_thomas2024_f)

cd8_thomas <- autoTune4(cd8_thomas)
cd8_thomas <- autoPilot(cd8_thomas, resolution = 0.5)
get_all_markers(cd8_thomas)

cd4_thomas <- subsetter(cd8_thomas, c(1, 3, 4, 5, 6, 7), original_obj_thomas2024_f)

cd4_thomas <- autoTune4(cd4_thomas)
cd4_thomas <- autoPilot(cd4_thomas, resolution = 0.5)
get_all_markers(cd4_thomas)

cd4.1_thomas <- subsetter(cd4_thomas, c(0, 1, 2, 5, 6, 7, 8, 9), original_obj_thomas2024_f)

cd4.1_thomas <- autoTune4(cd4.1_thomas)
cd4.1_thomas <- autoPilot(cd4.1_thomas, resolution = 0.5)
get_all_markers(cd4.1_thomas)


# 3 plotting
## 3.0 defining groups
seurat_obj <- cd4.1_thomas
cl <- as.character(seurat_obj$seurat_clusters)
cl[cl == "7"] <- "6"
seurat_obj$seurat_clusters <- factor(cl)
Idents(seurat_obj) <- seurat_obj$seurat_clusters
table(seurat_obj$seurat_clusters)

cluster_names <- c(
  "0" = "Resting",
  "1" = "Activated", 
  "2" = "Naive",
  "3" = "Treg",
  "4" = "IL-21 Effector",
  "5" = "Tfh",
  "6" = "Proliferating",
  "7" = ".",
  "8" = "Activated Treg"
)

seurat_obj@meta.data$cluster_names <- cluster_names[as.character(seurat_obj$seurat_clusters)]

# Define groups
inflamed <- c("SIC_100", "SIC_43", "SIC_71", "SIC_89", "SIC_97", "SIC_76", "SIC_121", "SIC_126", "SIC_134", "SIC_140", "SIC_141", "SIC_32", "SIC_36", "SIC_40")
control  <- c("SIC_13", "SIC_186", "SIC_187", "SIC_188", "SIC_196", "MC_1", "MC_2", "MC_9")
treated  <- c("SIC_109", "SIC_172", "SIC_19", "SIC_31", "SIC_53", "SIC_94", "SIC_33", "SIC_132")

seurat_obj@meta.data <- seurat_obj@meta.data %>%
  mutate(condi = case_when(
    sample %in% inflamed ~ 'IrAE',
    sample %in% control  ~ 'Control',
    sample %in% treated  ~ 'other',
    TRUE ~ 'unknown' 
  ))

## re-subsetting
obj <- seurat_obj

cells8 <- WhichCells(obj, idents = "8")
expr   <- FetchData(obj, vars = c("IL21", "FOXP3"), cells = cells8)
both   <- sum(expr$IL21 > 0 & expr$FOXP3 > 0)
cat("cluster 8 cells:", length(cells8),
    "| IL21+:", sum(expr$IL21 > 0),
    "| FOXP3+:", sum(expr$FOXP3 > 0),
    "| both:", both, "\n")
tab <- table(IL21 = expr$IL21 > 0, FOXP3 = expr$FOXP3 > 0)

## isolating cluster 8 (Tph/IL-21 effector)
cluster8 <- subset(seurat_obj, subset = seurat_clusters == "8")
cluster8 <- autoTune4(cluster8)
cluster8 <- autoPilot(cluster8, resolution = 0.1)

# Label each subcluster
cluster8$final_label <- dplyr::case_when(
  cluster8$seurat_clusters == "0" ~ "IL-21 Effector",
  cluster8$seurat_clusters == "1" ~ "Treg",
  TRUE                            ~ "Treg"
)

# Overwrite only the cluster-8 barcodes in the cluser name
new_labels <- setNames(cluster8$final_label, colnames(cluster8))
obj$cluster_names[names(new_labels)] <- new_labels

seurat_obj <- obj

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

## 3.2 Heatmap

key_markers <- c("PDCD1", "CXCR5", "CXCL13", "IL21", "IL23R", "IL17A", "RORC", "FOXP3", "TBX21", "MKI67", "GATA3", "IFNG", "CXCR3", "CCR6", "CCR7", "CD69")
avg_exp_scaled <- t(scale(t(AverageExpression(seurat_obj, features = key_markers, group.by = "seurat_clusters")$RNA)))

cluster_label_map <- setNames(
  seurat_obj@meta.data[["cluster_names"]],
  as.character(seurat_obj@meta.data[["seurat_clusters"]])
)
cluster_label_map <- cluster_label_map[!duplicated(names(cluster_label_map))]

# Strip "g" prefix from column names to match
cluster_label_map <- cluster_label_map[sub("^g", "", colnames(avg_exp_scaled))]

Heatmap(avg_exp_scaled,
        name = "Scaled\nExpression",
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        cluster_rows = FALSE, cluster_columns = FALSE,
        column_labels = cluster_label_map,
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 9))

## 3.3 UMAP
cluster_colors_cd4 <- c("0" = "#ADACCB", "1" = "#F9989A" , "2" =  "#68397A", 
                        "4" = "#F5CB1A", "3" = "#CB449F", "5" = "mediumseagreen",
                        "6" = "palevioletred", "7" = "skyblue4", "8" = "#F5CB1A",
                        "9" =  "skyblue2", "10" = "burlywood2", "11" = "#3D5EA9" , 
                        "12" = "khaki", "13" = "darkviolet", "14" = "turquoise", 
                        "15" = "#CBC8BB", "16" =  "royalblue"  , "17" = "#DF0000", 
                        "18" = "peru", "19" = "orange", "20" = "grey30",
                        "21" = "springgreen3" , "22" = "seagreen", "23" = "maroon4",
                        "24" = "black", "25" = "hotpink") 

plot_umap_labeled <- function(seurat_obj, label_col = "cluster_names", title = "CD4+ T Cell Clusters") {
  
  if (!label_col %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Column '%s' not found in meta.data. Available columns: %s",
                 label_col, paste(colnames(seurat_obj@meta.data), collapse = ", ")))
  }
  
  umap_df <- as.data.frame(Embeddings(seurat_obj, "umap"))
  colnames(umap_df) <- c("UMAP_1", "UMAP_2")
  umap_df$label <- seurat_obj@meta.data[[label_col]]
  umap_df$cluster_id <- as.character(seurat_obj@meta.data[["seurat_clusters"]])
  
  # Build color vector keyed by label using seurat_clusters as bridge
  label_to_color <- setNames(
    cluster_colors_cd4[umap_df$cluster_id],
    umap_df$label
  )
  label_to_color <- label_to_color[!duplicated(names(label_to_color))]
  
  centroids <- umap_df %>%
    group_by(label) %>%
    summarise(across(starts_with("UMAP"), median))
  
  ggplot(umap_df, aes(UMAP_1, UMAP_2, color = label)) +
    geom_point(size = 0.5, alpha = 0.9, stroke = 0) +
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
seurat_obj$origin <- "thomas2024"
seurat_obj$pt <- seurat_obj$sample
pt_col        <- "pt"                   
condition_col <- "condi"            
dataset_col   <- "origin" 
cluster_col   <- "cluster_names"     

conditions_to_compare <- c("Control", "IrAE")

clusters_to_plot <- NULL

cond_fill  <- c("#F8C5C8", "#E63946")  
cond_point <- c("#E89AA0", "#C81E2C")   

names(cond_fill)  <- conditions_to_compare
names(cond_point) <- conditions_to_compare

md <- seurat_obj@meta.data %>%
  dplyr::select(all_of(c(pt_col, condition_col, dataset_col, cluster_col))) %>%
  dplyr::rename(pt = !!pt_col,
                condition = !!condition_col,
                dataset = !!dataset_col,
                cluster = !!cluster_col) %>%
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

plot_cluster <- function(cl) {
  df <- filter(plot_df, cluster == cl)
  ggplot(df, aes(x = condition, y = pct, fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.6) +
    geom_jitter(aes(color = condition), width = 0.15, size = 2, alpha = 0.85) +
    scale_fill_manual(values = cond_fill) +
    scale_color_manual(values = cond_point) +
    scale_y_continuous(limits = c(0, 50), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Condition",
         y = "% of total cells (per sample)",
         title = paste("Cluster", cl)) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none")
}

cls <- if (is.null(clusters_to_plot)) sort(unique(plot_df$cluster)) else clusters_to_plot
plots <- lapply(cls, plot_cluster)
names(plots) <- cls

plots

## 3.5 statistics
sample_col    <- "sample"
condition_col <- "condi"
cluster_col   <- "cluster_names"

cond_A        <- "control"   
cond_B        <- "inflamed"   

meta <- seurat_obj@meta.data

meta <- meta[meta[[condition_col]] %in% c(cond_A, cond_B), ]

prop_table <- meta %>%
  dplyr::count(.sample = .data[[sample_col]],
               .condition = .data[[condition_col]],
               .cluster = .data[[cluster_col]]) %>%
  group_by(.sample) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

prop_table %>%
  group_by(.cluster) %>%
  summarise(p_value = wilcox.test(proportion ~ .condition)$p.value) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))
















