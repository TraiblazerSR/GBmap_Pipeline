# =============================================================================
# Step 6: Pseudotime Analysis with Monocle3
# =============================================================================
# This script performs pseudotime analysis using monocle3 for target cell types
# and target gene co-trajectory analysis. 
# =============================================================================

# Clear workspace and set working directory
rm(list = ls())
project_name = "GBmap_Pipeline_v1.0"
if (getwd() != "/data/y1005/GBmap_Pipeline_v1.0") {
  setwd("/data/y1005/GBmap_Pipeline_v1.0")
}

# Load global config and seed
source("config/common_utils.R")
cfg <- read_global_config()
set.seed(cfg$seed)

options(future.globals.maxSize = 50 * 1024^3)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(Matrix)
  library(monocle3)
  library(future)
  library(pheatmap)
  library(viridis)
})

future::plan("multisession", workers = cfg$workers)

# =============================================================================
# SECTION 1: CONFIGURATION & PATHS ####

# Load target cell types from config (step6 specific)
target_cell_types <- if (!is.null(cfg$step6_target_cell_types)) {
  cfg$step6_target_cell_types
} else {
  cfg$target_cell_types
}

# Input/output paths
input_rds_path <- "step1_add_annotations_from_crosstable/processed/core_data.rds"
if (!file.exists(input_rds_path)) stop("Input RDS not found: ", input_rds_path)

cross_table_path <- cfg$annotation_cross_table_monocle_path
if (!file.exists(cross_table_path)) stop("Cross table not found: ", cross_table_path)

step6_root_dir <- "step6_pseudotime_target_gene"
results_dir <- file.path(step6_root_dir, "results")
figures_dir <- file.path(step6_root_dir, "figures")
processed_dir <- file.path(step6_root_dir, "processed")

for (d in c(results_dir, figures_dir, processed_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# =============================================================================
# SECTION 2: UTILITY FUNCTIONS ####

save_gg <- function(plot_obj, file_base, width = 7, height = 5) {
  ggplot2::ggsave(filename = paste0(file_base, ".pdf"), plot = plot_obj, width = width, height = height, limitsize = FALSE)
  ggplot2::ggsave(filename = paste0(file_base, ".png"), plot = plot_obj, width = width, height = height, dpi = 300, limitsize = FALSE)
}

get_counts_matrix <- function(seurat_object, assay_name = "RNA") {
  counts_mat <- NULL
  try({
    counts_mat <- SeuratObject::GetAssayData(seurat_object, assay = assay_name, layer = "counts")
  }, silent = TRUE)
  if (is.null(counts_mat) || length(counts_mat) == 0) {
    counts_mat <- SeuratObject::GetAssayData(seurat_object, assay = assay_name, slot = "counts")
  }
  if (!inherits(counts_mat, "dgCMatrix")) counts_mat <- as(as.matrix(counts_mat), "dgCMatrix")
  counts_mat
}

join_annotation_by_level4 <- function(seurat_object, cross_table, out_col = "annotation_monocle") {
  if (!("annotation_level_4" %in% colnames(seurat_object@meta.data))) {
    stop("annotation_level_4 not found in Seurat meta.data. Ensure Step1 data is used.")
  }
  meta <- seurat_object@meta.data %>% tibble::rownames_to_column(".cell")
  key <- cross_table %>% dplyr::select(annotation_level_4, annotation_custom) %>% dplyr::distinct()
  out <- meta %>% dplyr::left_join(key, by = "annotation_level_4")
  if (!("annotation_custom" %in% colnames(out))) stop("annotation_custom column missing after join.")
  seurat_object@meta.data[[out_col]] <- out$annotation_custom
  seurat_object
}

resolve_target_feature <- function(seurat_object, primary_symbol, aliases_csv = NULL, ensembl_id = NULL) {
  candidates <- unique(na.omit(c(primary_symbol,
                                 if (!is.null(aliases_csv)) strsplit(aliases_csv, ",")[[1]] else NULL,
                                 ensembl_id)))
  candidates <- trimws(candidates)
  features <- rownames(seurat_object)
  # 1) exact match to feature rownames
  hit <- intersect(candidates, features)
  if (length(hit) > 0) return(hit[[1]])
  # 2) try to locate by meta feature names if present
  mf <- try(seurat_object@assays$RNA@meta.features, silent = TRUE)
  if (!inherits(mf, "try-error") && nrow(mf) > 0) {
    for (col in c("gene_short_name", "gene_name", "feature_name")) {
      if (col %in% colnames(mf)) {
        idx <- which(!is.na(mf[[col]]) & mf[[col]] %in% candidates)
        if (length(idx) > 0) return(rownames(mf)[idx[[1]]])
      }
    }
  }
  stop("Target gene not found in Seurat features: ", paste(candidates, collapse = ", "))
}

build_cds_from_seurat <- function(seurat_object, assay_name = "RNA") {
  counts_mat <- get_counts_matrix(seurat_object, assay_name = assay_name)
  cell_metadata <- seurat_object@meta.data
  gene_metadata <- data.frame(gene_short_name = rownames(counts_mat), row.names = rownames(counts_mat), stringsAsFactors = FALSE)
  monocle3::new_cell_data_set(expression_data = counts_mat, cell_metadata = cell_metadata, gene_metadata = gene_metadata)
}

preprocess_and_learn <- function(cds) {
  cds <- monocle3::preprocess_cds(cds, num_dim = 100)
  cds <- monocle3::reduce_dimension(cds)
  cds <- monocle3::cluster_cells(cds)
  cds <- monocle3::learn_graph(cds, use_partition = FALSE, close_loop = FALSE)
  cds
}

compute_correlations_against_target <- function(cds, target_rowname) {
  pt <- monocle3::pseudotime(cds)
  finite_cells <- is.finite(pt)
  if (!any(finite_cells)) stop("No finite pseudotime cells found after ordering.")
  expr <- monocle3::exprs(cds)[, finite_cells]
  target_vec <- as.numeric(monocle3::exprs(cds)[target_rowname, finite_cells])
  cors <- apply(expr, 1, function(g) suppressWarnings(stats::cor(g, target_vec, method = "spearman", use = "complete.obs")))
  
  # Handle rowData structure
  row_data <- SummarizedExperiment::rowData(cds)
  if (is.data.frame(row_data) && "gene_short_name" %in% colnames(row_data)) {
    gene_short_names <- row_data[names(cors), ]$gene_short_name
  } else {
    gene_short_names <- if (is.atomic(row_data)) {
      row_data[names(cors)]
    } else {
      names(cors)
    }
  }
  
  tibble::tibble(
    gene_id = names(cors),
    gene_short_name = gene_short_names,
    spearman_rho = as.numeric(cors)
  ) %>% dplyr::arrange(desc(abs(spearman_rho)))
}

run_graph_test_and_modules <- function(cds) {
  res <- monocle3::graph_test(cds, neighbor_graph = "principal_graph", cores = cfg$workers)
  res$gene_id <- rownames(res)
  sig_ids <- rownames(subset(res, q_value < 0.05))
  if (length(sig_ids) == 0) {
    return(list(graph_test = res, gene_module_df = NULL))
  }
  gene_module_df <- monocle3::find_gene_modules(cds[sig_ids, ], resolution = 1e-2, cores = cfg$workers)
  return(list(graph_test = res, gene_module_df = gene_module_df))
}

# =============================================================================
# SECTION 3: LOAD DATA & ANNOTATIONS ####

cat("Loading Seurat object...\n")
seurat_obj <- readRDS(input_rds_path)
cat("Loaded Seurat object with", ncol(seurat_obj), "cells and", nrow(seurat_obj), "features\n")

cat("Loading annotation cross-table...\n")
cross_table <- readr::read_csv(cross_table_path, show_col_types = FALSE)
cat("Cross-table dimensions:", paste(dim(cross_table), collapse = " x "), "\n")

# Add monocle annotation
seurat_obj <- join_annotation_by_level4(seurat_obj, cross_table, out_col = "annotation_monocle")

# Resolve target gene
target_gene_rowname <- resolve_target_feature(
  seurat_obj, 
  primary_symbol = cfg$target_gene_symbol,
  aliases_csv = cfg$target_gene_aliases,
  ensembl_id = cfg$target_gene_ensembl
)
cat("Target gene resolved to:", target_gene_rowname, "\n")

# =============================================================================
# SECTION 4: PSEUDOTIME ANALYSIS PER CELL TYPE ####

for (cell_type in target_cell_types) {
  cat("\n=== Processing cell type:", cell_type, "===\n")
  
  # Create cell type specific directories
  ct_results_dir <- file.path(results_dir, gsub("[^A-Za-z0-9]+", "_", cell_type))
  ct_figures_dir <- file.path(figures_dir, gsub("[^A-Za-z0-9]+", "_", cell_type))
  
  for (d in c(ct_results_dir, ct_figures_dir)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  
  # Filter cells for this cell type
  cell_mask <- seurat_obj@meta.data$annotation_monocle == cell_type
  if (sum(cell_mask, na.rm = TRUE) < 10) {
    cat("Skipping", cell_type, "- insufficient cells (", sum(cell_mask, na.rm = TRUE), ")\n")
    next
  }
  
  seurat_ct <- seurat_obj[, cell_mask]
  cat("Cell type subset:", cell_type, "cells =", ncol(seurat_ct), "\n")
  
  # Build monocle3 CDS
  cat("Building monocle3 CDS...\n")
  cds <- build_cds_from_seurat(seurat_ct)
  
  # Preprocess and learn trajectory
  cat("Preprocessing and learning trajectory...\n")
  cds <- preprocess_and_learn(cds)
  
  # Order cells (automatic root selection - lowest target gene expression)
  cat("Ordering cells along pseudotime...\n")
  target_expr <- monocle3::exprs(cds)[target_gene_rowname, ]
  root_cells <- names(target_expr)[which.min(target_expr)]
  cds <- monocle3::order_cells(cds, root_cells = root_cells)
  
  # Save CDS
  saveRDS(cds, file.path(ct_results_dir, paste0("cds_", gsub("[^A-Za-z0-9]+", "_", cell_type), ".rds")))
  
  # Generate plots
  cat("Generating pseudotime plots...\n")
  
  # UMAP with pseudotime
  p_umap_pt <- monocle3::plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                                   cell_size = 0.5) +
    scale_color_viridis_c(name = "Pseudotime") +
    ggtitle(paste("Pseudotime -", cell_type)) +
    theme_minimal()
  save_gg(p_umap_pt, file.path(ct_figures_dir, paste0("pseudotime_umap_", gsub("[^A-Za-z0-9]+", "_", cell_type))))
  
  # UMAP with target gene expression
  p_umap_target <- monocle3::plot_cells(cds, genes = target_gene_rowname, label_cell_groups = FALSE,
                                       cell_size = 0.5) +
    scale_color_viridis_c(name = paste(cfg$target_gene_symbol, "Expression")) +
    ggtitle(paste(cfg$target_gene_symbol, "Expression -", cell_type)) +
    theme_minimal()
  save_gg(p_umap_target, file.path(ct_figures_dir, paste0("target_gene_umap_", gsub("[^A-Za-z0-9]+", "_", cell_type))))
  
  # Trajectory plot
  p_trajectory <- monocle3::plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
                                      show_trajectory_graph = TRUE, cell_size = 0.5) +
    scale_color_viridis_c(name = "Pseudotime") +
    ggtitle(paste("Trajectory -", cell_type)) +
    theme_minimal()
  save_gg(p_trajectory, file.path(ct_figures_dir, paste0("trajectory_", gsub("[^A-Za-z0-9]+", "_", cell_type))))
  
  # Compute correlations with target gene
  cat("Computing correlations with target gene...\n")
  correlations <- compute_correlations_against_target(cds, target_gene_rowname)
  readr::write_csv(correlations, file.path(ct_results_dir, paste0("target_gene_correlations_", gsub("[^A-Za-z0-9]+", "_", cell_type), ".csv")))
  
  # Top correlated genes heatmap
  top_genes <- correlations %>% 
    dplyr::filter(!is.na(spearman_rho)) %>%
    dplyr::slice_head(n = 50)
  
  if (nrow(top_genes) > 5) {
    cat("Creating correlation heatmap...\n")
    pt <- monocle3::pseudotime(cds)
    finite_cells <- is.finite(pt)
    expr_mat <- monocle3::exprs(cds)[top_genes$gene_id, finite_cells]
    
    # Order cells by pseudotime
    pt_ordered <- sort(pt[finite_cells])
    expr_mat_ordered <- expr_mat[, names(pt_ordered)]
    
    # Create heatmap
    pdf(file.path(ct_figures_dir, paste0("correlation_heatmap_", gsub("[^A-Za-z0-9]+", "_", cell_type), ".pdf")), width = 12, height = 8)
    pheatmap::pheatmap(expr_mat_ordered, 
                      cluster_cols = FALSE, 
                      cluster_rows = TRUE,
                      scale = "row",
                      color = viridis::viridis(100),
                      show_colnames = FALSE,
                      main = paste("Top Target-Correlated Genes -", cell_type))
    dev.off()
    
    png(file.path(ct_figures_dir, paste0("correlation_heatmap_", gsub("[^A-Za-z0-9]+", "_", cell_type), ".png")), width = 12, height = 8, units = "in", res = 300)
    pheatmap::pheatmap(expr_mat_ordered, 
                      cluster_cols = FALSE, 
                      cluster_rows = TRUE,
                      scale = "row",
                      color = viridis::viridis(100),
                      show_colnames = FALSE,
                      main = paste("Top Target-Correlated Genes -", cell_type))
    dev.off()
  }
  
  # Graph test and gene modules
  cat("Running graph test and finding gene modules...\n")
  graph_results <- run_graph_test_and_modules(cds)
  
  # Save results
  readr::write_csv(graph_results$graph_test, file.path(ct_results_dir, paste0("graph_test_", gsub("[^A-Za-z0-9]+", "_", cell_type), ".csv")))
  
  if (!is.null(graph_results$gene_module_df)) {
    readr::write_csv(graph_results$gene_module_df, file.path(ct_results_dir, paste0("gene_modules_", gsub("[^A-Za-z0-9]+", "_", cell_type), ".csv")))
  }
  
  cat("Completed analysis for", cell_type, "\n")
}

# =============================================================================
# SECTION 5: SESSION INFO & SUMMARY ####

session_info <- sessionInfo()
saveRDS(session_info, file.path(results_dir, "session_info.rds"))

# Create summary
analysis_summary <- data.frame(
  Cell_Type = target_cell_types,
  Cells_Analyzed = sapply(target_cell_types, function(ct) {
    mask <- seurat_obj@meta.data$annotation_monocle == ct
    sum(mask, na.rm = TRUE)
  }),
  Target_Gene = cfg$target_gene_symbol,
  Target_Gene_ID = target_gene_rowname,
  Analysis_Date = Sys.Date()
)

readr::write_csv(analysis_summary, file.path(results_dir, "analysis_summary.csv"))

cat("\n=== Analysis Complete ===\n")
cat("Results saved in:", results_dir, "\n")
cat("Figures saved in:", figures_dir, "\n")
cat("Summary:\n")
print(analysis_summary)

cat("\nScript completed successfully!\n")
