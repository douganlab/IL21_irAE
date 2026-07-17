# EC_functions.R
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

# Enoch's Functions

subsetter <- function(seurat_obj, clusters, original_obj) { 
  clusters_to_extract <- clusters
  cells_to_keep <- WhichCells(seurat_obj, idents = clusters_to_extract)
  
  unprocessed_subset <- subset(original_obj, cells = cells_to_keep)
}

cd8_subsetter <- function(seurat_obj, original_obj, threshold = 0) { 
  cd8_cells <- WhichCells(seurat_obj, expression = CD8A > threshold | CD8B > threshold)
  all_cells <- colnames(seurat_obj)
  cells_to_keep <- setdiff(all_cells, cd8_cells)
  unprocessed_subset <- subset(original_obj, cells = cells_to_keep)
  return(unprocessed_subset)
}

merge_by_common_genes <- function(seurat_list, output_name) {
  gene_names_list <- lapply(seurat_list, function(x) rownames(x))
  cat("Number of genes in each dataset:\n")
  for (i in 1:length(gene_names_list)) {
    cat("Dataset", i, ":", length(gene_names_list[[i]]), "genes\n")
  }
  common_genes <- gene_names_list[[1]]
  for (i in 2:length(gene_names_list)) {
    common_genes <- intersect(common_genes, gene_names_list[[i]])
  }
  cat("\nNumber of common genes across all datasets:", length(common_genes), "\n")
  merged_obj <- merge(seurat_list[[1]][common_genes, ], 
                      y = lapply(seurat_list[-1], function(x) x[common_genes, ]))
  new_name <- paste0("original_obj_", output_name)
  assign(new_name, merged_obj, envir = .GlobalEnv)
  assign(output_name, merged_obj, envir = .GlobalEnv)
  cat("\nCreated objects:", output_name, "and", new_name, "\n")
  return(invisible(merged_obj))
}

# Merge by union ---- 
merge_union <- function(seurat_list, output_name) {
  cat("Merging", length(seurat_list), "samples by union...\n")
  all_genes <- unique(unlist(lapply(seurat_list, rownames)))
  cat("Number of genes in each dataset:\n")
  for (i in 1:length(seurat_list)) {
    cat("Dataset", i, ":", nrow(seurat_list[[i]]), "genes\n")
  }
  cat("\nNumber of total unique genes:", length(all_genes), "\n")
  expand_genes <- function(seurat_obj, all_genes) {
    current_genes <- rownames(seurat_obj)
    missing_genes <- setdiff(all_genes, current_genes)
    
    if (length(missing_genes) > 0) {
      counts <- GetAssayData(seurat_obj, layer = "counts")
      
      zero_matrix <- Matrix::Matrix(0, 
                                    nrow = length(missing_genes), 
                                    ncol = ncol(counts),
                                    sparse = TRUE)
      rownames(zero_matrix) <- missing_genes
      colnames(zero_matrix) <- colnames(counts)
      
      expanded_counts <- rbind(counts, zero_matrix)
      
      expanded_counts <- expanded_counts[all_genes, ]
      
      seurat_obj_new <- CreateSeuratObject(
        counts = expanded_counts,
        meta.data = seurat_obj@meta.data
      )
      
      return(seurat_obj_new)
    }
    
    return(seurat_obj)
  }
  cat("Expanding all samples to include all genes...\n")
  expanded_list <- lapply(seurat_list, expand_genes, all_genes = all_genes)
  cat("Merging expanded samples...\n")
  merged_obj <- merge(
    x = expanded_list[[1]],
    y = expanded_list[2:length(expanded_list)],
    merge.data = FALSE
  )
  cat("Merge complete:", ncol(merged_obj), "cells,", nrow(merged_obj), "genes\n")
  new_name <- paste0("original_obj_", output_name)
  assign(new_name, merged_obj, envir = .GlobalEnv)
  assign(output_name, merged_obj, envir = .GlobalEnv)
  cat("\nCreated objects:", output_name, "and", new_name, "\n")
  return(invisible(merged_obj))
}

merge_union2 <- function(seurat_list, output_name, cell_ids = NULL) {
  cat("Merging", length(seurat_list), "samples by union...\n")
  all_genes <- unique(unlist(lapply(seurat_list, rownames)))
  cat("Number of genes in each dataset:\n")
  for (i in 1:length(seurat_list)) {
    cat("Dataset", i, ":", nrow(seurat_list[[i]]), "genes\n")
  }
  cat("\nNumber of total unique genes:", length(all_genes), "\n")
  expand_genes <- function(seurat_obj, all_genes) {
    seurat_obj <- JoinLayers(seurat_obj)
    current_genes <- rownames(seurat_obj)
    missing_genes <- setdiff(all_genes, current_genes)
    if (length(missing_genes) > 0) {
      counts <- GetAssayData(seurat_obj, layer = "counts")
      
      zero_matrix <- Matrix::Matrix(0,
                                    nrow = length(missing_genes),
                                    ncol = ncol(counts),
                                    sparse = TRUE)
      rownames(zero_matrix) <- missing_genes
      colnames(zero_matrix) <- colnames(counts)
      
      expanded_counts <- rbind(counts, zero_matrix)[all_genes, ]
      
      seurat_obj <- CreateSeuratObject(
        counts = expanded_counts,
        meta.data = seurat_obj@meta.data
      )
    }
    return(seurat_obj)
  }
  expanded_list <- lapply(seurat_list, expand_genes, all_genes = all_genes)
  ids <- if (!is.null(cell_ids)) cell_ids else if (!is.null(names(seurat_list))) names(seurat_list) else NULL
  merged <- merge(
    x = expanded_list[[1]],
    y = expanded_list[-1],
    add.cell.ids = ids,
    project = output_name
  )
  merged <- JoinLayers(merged)
  cat("\nFinal merged object:", ncol(merged), "cells x", nrow(merged), "genes\n")
  return(merged)
}

autoTune4 <- function(seurat_obj, 
                      min_features = 200, 
                      max_features = 6000, 
                      max_mt = 25,
                      n_variable_features = 3000,
                      n_pcs = 50,
                      harmony_dims = 1:30,
                      batch_var = "sample") {
  seurat_obj <- JoinLayers(seurat_obj)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  cat("Step 1: Mitochondrial exclusion and filtering...\n")
  seurat_obj <- subset(seurat_obj, 
                       subset = nFeature_RNA > min_features & 
                         nFeature_RNA < max_features & 
                         percent.mt < max_mt)
  cat("  Cells remaining after filtering:", ncol(seurat_obj), "\n\n")
  cat(ncol(seurat_obj), " cells \n")
  cat(nrow(seurat_obj), " genes \n")
  filtered_genes <- rownames(seurat_obj)
  filtered_genes <- filtered_genes[!grepl("^MT-", filtered_genes)]
  filtered_genes <- filtered_genes[rowSums(LayerData(seurat_obj, layer = "counts")[filtered_genes, ]) > 0]
  seurat_obj <- seurat_obj[filtered_genes, ]
  cat(ncol(seurat_obj), " cells \n")
  cat(nrow(seurat_obj), " genes \n")
  cat("Step 2: Normalize and find variable genes...\n")
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, 
                                     selection.method = "vst", 
                                     nfeatures = n_variable_features)
  cat("  Found", length(VariableFeatures(seurat_obj)), "variable features\n\n")
  cat("Step 3: Scale data and run PCA...\n")
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, npcs = n_pcs, verbose = FALSE)
  cat("  PCA complete with", n_pcs, "PCs\n\n")
  cat("Step 4: Run Harmony batch correction...\n")
  seurat_obj <- RunHarmony(
    object = seurat_obj,
    group.by.vars = batch_var,
    dims = harmony_dims,
    verbose = TRUE
  )
  cat("\nAutoTune complete!\n")
  return(seurat_obj)
}


autoPilot <- function(seurat_obj, 
                      resolution,
                      harmony_dims = 1:30,
                      reduction = "harmony", 
                      k.param = 20) {
  
  cat("Step 1: Running UMAP...\n")
  seurat_obj <- RunUMAP(seurat_obj, reduction = reduction, dims = harmony_dims)
  cat("  UMAP complete\n\n")
  cat("Step 2: Finding neighbors...\n")
  seurat_obj <- FindNeighbors(seurat_obj, reduction = reduction, dims = harmony_dims, k.param = k.param)
  cat("  Neighbors found\n\n")
  cat("Step 3: Finding clusters at resolution", resolution, "...\n")
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  cat("  Clustering complete\n\n")
  cat("autoMapper complete!\n")
  return(seurat_obj)
}

get_top50_markers_log <- function(seurat_obj) {
  obj_name <- deparse(substitute(seurat_obj))
  allMarkers <- FindAllMarkers(seurat_obj, only.pos = TRUE, logfc.threshold = 0.25)
  top50_per_cluster <- allMarkers %>%
    group_by(cluster) %>%
    top_n(n = 50, wt = avg_log2FC) %>%
    ungroup()
  filename <- paste0(obj_name, "_markers.csv")
  write.csv(top50_per_cluster, filename, row.names = FALSE)
  cat("Saved markers to:", filename, "\n")
  return(top50_per_cluster)
}


get_top50_markers <- function(seurat_obj) {
  markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    test.use = "wilcox"
  )
  
  markers %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = 50) %>%
    ungroup() %>%
    arrange(cluster, desc(avg_log2FC))
}


get_all_markers <- function(seurat_obj) {
  obj_name <- deparse(substitute(seurat_obj))
  allMarkers <- FindAllMarkers(seurat_obj, only.pos = TRUE, logfc.threshold = 0.25)
  filename <- paste0(obj_name, "_markers.csv")
  write.csv(allMarkers, filename, row.names = FALSE)
  cat("Saved markers to:", filename, "\n")
  return(allMarkers)
}
