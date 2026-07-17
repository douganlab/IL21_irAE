# combined.R
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
## Luoma
luoma <- readRDS('/cd4.1_luoma.rds')
luoma_o <- readRDS('/original_obj_luoma2020_f.rds')
## Thomas 
thomas <- readRDS('/cd4.1_thomas.rds')
thomas_o <- readRDS('/original_obj_thomas2024_f.rds')
## Lechner
lechner <- readRDS('/cd4.3_lechner.rds')
lechner_o <- readRDS('/original_obj_lechner2023_f.rds')
## Zhang
zhang <- readRDS('/seurat_obj.rds')
zhang_o <- readRDS('/original_obj_zhang2024_f.rds')
## Blum
blum <- readRDS('/seurat_obj_blum.rds')
blum_o <- readRDS('/original_obj_blum2024_f.rds')

## 1.4 Metatdata addition
# metadata coded per the publically available data 
luoma$condi <- case_when(
  luoma$pt %in% c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8") ~ "IrAE",
  luoma$pt %in% c("CT1", "CT2", "CT3", "CT4", "CT5", "CT6", "CT7", "CT8") ~ "Control",
  luoma$pt %in% c("NC1", "NC2", "NC3", "NC4", "NC5", "NC6") ~ "no_IrAE",
  TRUE ~ NA_character_
)
thomas$condi <- case_when(
  thomas$pt %in% c("SIC_100", "SIC_43", "SIC_71", "SIC_89", "SIC_97", "SIC_76", "SIC_121", "SIC_126", "SIC_134", "SIC_140", "SIC_141", "SIC_32", "SIC_36", "SIC_40") ~ "IrAE",
  thomas$pt %in% c("SIC_13", "SIC_186", "SIC_187", "SIC_188", "SIC_196", "MC_1", "MC_2", "MC_9") ~ "Control",
  thomas$pt %in% c("SIC_109", "SIC_172", "SIC_19", "SIC_31", "SIC_53", "SIC_94", "SIC_33", "SIC_132") ~ "no_IrAE",
  TRUE ~ NA_character_
)
lechner$condi <- case_when(
  lechner$condition == "IRAE" ~ "IrAE",
  lechner$condition == "Normal" ~ "Control",
  lechner$condition == "HT" ~ "Other",
  TRUE ~ NA_character_
)
zhang$condi <- case_when(
  zhang$pt %in% c("BALF-1-AE", "BALF-2-AE", "BALF-3-AE", "BALF-4-AE", "BALF-5-AE") ~ "IrAE",
  zhang$pt %in% c("C51", "C52", "C100") ~ "Control",
  TRUE ~ NA_character_
)
blum$condi <- case_when(
  blum$pt %in% c("264", "258", "199", "48", "197", "177", "232", "164", "237", "217", "153", "175", "171", "317", "319") ~ "IrAE",
  blum$pt %in% c("182", "333") ~ "no_IrAE",
  TRUE ~ NA_character_
)
# saved here
saveRDS(luoma, "luoma.rds")
saveRDS(thomas, "thomas.rds")
saveRDS(lechner, "lechner.rds")
saveRDS(zhang, "zhang.rds")
saveRDS(blum, "blum.rds")
# we repeated above with the "_o" variants and also saved them. 
saveRDS(luoma_o, "luoma_o.rds")
saveRDS(thomas_o, "thomas_o.rds")
saveRDS(lechner_o, "lechner_o.rds")
saveRDS(zhang_o, "zhang_o.rds")
saveRDS(blum_o, "blum_o.rds")

## 1.5 Merging
combined <- merge_union2(seurat_list = list(luoma = luoma, thomas = thomas, lechner = lechner, zhang = zhang, blum = blum), output_name = "combined")
original_combined <- merge_union2(seurat_list = list(luoma = luoma_o, thomas = thomas_o, lechner = lechner_o, zhang = zhang_o, blum = blum_o), output_name = "original_combined")

# 2 Processing
combined <- autoTune4(combined)
combined <- autoPilot(combined, resolution = 1)
get_all_markers(combined)

cd4.0 <- subsetter(combined, c(0, 2, 3, 4, 5, 6, 7, 8, 10, 11, 13, 14, 15, 16, 18), original_combined)

cd4.0 <- autoTune4(cd4.0)
cd4.0 <- autoPilot(cd4.0, resolution = 0.5)
get_all_markers(cd4.0)

cd4.1 <- subsetter(cd4.0, c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), original_combined)

cd4.1 <- autoTune4(cd4.1)
cd4.1 <- autoPilot(cd4.1, resolution = 0.5)
get_all_markers(cd4.1)

# 3 Plotting
## UMAP
seurat_obj <- cd4.1
cl <- as.character(seurat_obj$seurat_clusters)
cl[cl == "9"] <- "7"
cl[cl == "6"] <- "2"
seurat_obj$seurat_clusters <- factor(cl)
Idents(seurat_obj) <- seurat_obj$seurat_clusters
table(seurat_obj$seurat_clusters)

cluster_names <- c(
  "0" = "Activated",
  "1" = "Resting", 
  "2" = "Memory",
  "3" = "Treg",
  "4" = "Tfh", 
  "5" = "Th17",
  "6" = ".", 
  "7" = "Proliferating",
  "8" = "Tph",
  "9" = "..",
  "10" = "Exhausted",
  "11" = "Effector Treg"
)
seurat_obj@meta.data$cluster_names <- cluster_names[as.character(seurat_obj$seurat_clusters)]

seurat_obj$organ <- case_when(
  seurat_obj$origin %in% c("luoma2020", "thomas2024") ~ "Colon",
  seurat_obj$origin %in% "lechner2023" ~ "Thyroid",
  seurat_obj$origin %in% "zhang2024" ~ "Lung",
  seurat_obj$origin %in% "blum2024" ~ "heart",
  TRUE ~ NA_character_
)

cluster_colors_cd4 <- c("0" = "#ADACCB", "1" = "#F9989A" , "2" =  "#68397A", 
                        "4" = "#F5CB1A", "3" = "#CB449F", "5" = "mediumseagreen",
                        "6" = "palevioletred", "7" = "skyblue4", "8" = "brown",
                        "9" =  "skyblue2", "10" = "burlywood2", "11" = "#3D5EA9" , 
                        "12" = "khaki", "13" = "darkviolet", "14" = "turquoise", 
                        "15" = "#CBC8BB", "16" =  "royalblue"  , "17" = "#DF0000", 
                        "18" = "peru", "19" = "orange", "20" = "grey30",
                        "21" = "springgreen3" , "22" = "seagreen", "23" = "maroon4",
                        "24" = "black", "25" = "hotpink") 

## IL-21 plot
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

## Heatmap
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

## UMAP
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
    geom_point(size = 0.5, alpha = 0.6, stroke = 0) +
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

## expansion plot 
pt_col        <- "pt"            
condition_col <- "condi"       
dataset_col   <- "organ"             
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

stopifnot(all(conditions_to_compare %in% unique(as.character(md$condition))))

totals <- md %>%
  group_by(dataset, pt, condition) %>%
  summarise(total_cells = n(), .groups = "drop")

counts <- md %>%
  group_by(dataset, pt, condition, cluster) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  complete(cluster, nesting(dataset, pt, condition), fill = list(n_cells = 0))

plot_df <- counts %>%
  left_join(totals, by = c("dataset", "pt", "condition")) %>%
  mutate(pct = 100 * n_cells / total_cells,
         condition = factor(condition, levels = conditions_to_compare))

plot_cluster <- function(cl) {
  df <- filter(plot_df, cluster == cl)
  ggplot(df, aes(x = condition, y = pct, fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.6) +
    geom_jitter(aes(color = condition), width = 0.15, size = 2, alpha = 0.85) +
    facet_wrap(~ dataset, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = cond_fill) +
    scale_color_manual(values = cond_point) +
    scale_y_continuous(limits = c(0, 50), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Condition",
         y = "% of total cells (per organ)",
         title = paste("Cluster", cl)) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"))
}

cls <- if (is.null(clusters_to_plot)) sort(unique(plot_df$cluster)) else clusters_to_plot
plots <- lapply(cls, plot_cluster)
names(plots) <- cls

plots






























