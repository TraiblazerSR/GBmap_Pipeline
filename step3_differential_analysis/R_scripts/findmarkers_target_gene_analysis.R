# Clear workspace and set working directory
rm(list = ls())

# =============================================================================
# SECTION 1: SETUP & LOAD CONFIG ####
# =============================================================================
project_name <- "GBmap_Pipeline_v1.0"
if (getwd() != "/data/y1005/GBmap_Pipeline_v1.0") {
  setwd("/data/y1005/GBmap_Pipeline_v1.0")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(purrr)
})

source("config/common_utils.R")
cfg <- read_global_config()
set.seed(cfg$seed)

# Output directories (per-method will be created later)
ensure_directories(c(
  "step3_differential_analysis/processed",
  "step3_differential_analysis/results",
  "step3_differential_analysis/figures",
  "logs"
))

# Determine method order
method_setting <- tolower(cfg$group_methods)
method_order <- switch(method_setting,
  gmm = c("GMM"),
  median = c("median"),
  both = if (toupper(cfg$both_primary_method) == "GMM") c("GMM", "median") else c("median", "GMM"),
  c("GMM")
)

# =============================================================================
# SECTION 2: LOAD DATA AND GENE MAPPING ####
# =============================================================================

cat("Loading processed data from Step1...\n")
core_data <- readRDS("step1_add_annotations_from_crosstable/processed/core_data.rds")

# Use annotation_level_5 explicitly as in reference
if (!"annotation_level_5" %in% colnames(core_data@meta.data)) {
  stop("annotation_level_5 not found in core_data@meta.data")
}
all_cell_types_level <- "annotation_level_5"
all_cell_types <- unique(core_data@meta.data[[all_cell_types_level]])

# Identify target feature robustly
feature_name <- identify_target_gene_feature(core_data, cfg)
cat("Using target feature:", feature_name, "for symbol:", cfg$target_gene_symbol, "\n")

# Map genes for symbol annotations (ENSEMBL -> SYMBOL)
all_genes <- rownames(core_data)
tryCatch({
  gene_mapping <- clusterProfiler::bitr(all_genes, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  colnames(gene_mapping) <- c("ensembl_id", "gene_symbol")
}, error = function(e) {
  gene_mapping <<- data.frame(ensembl_id = character(0), gene_symbol = character(0))
})

# =============================================================================
# SECTION 3: PARAMETERS ####
# =============================================================================

analysis_method <- "FindMarkers"
min_cells_per_type <- 50
min_cells_per_group <- 10
logfc_threshold <- 0.25
min_pct <- 0.1
padj_threshold <- 0.05

primary_cell_types <- cfg$target_cell_types

# =============================================================================
# SECTION 4: CLASSIFICATION HELPERS ####
# =============================================================================

classify_cells <- function(expr_vec, method_label, cell_type_name) {
  if (method_label == "GMM") {
    res <- split_by_gmm(expr_vec)
    high <- res$high; low <- res$low
  } else {
    res <- split_by_median(expr_vec)
    high <- res$high; low <- res$low
  }
  list(high = high, low = low, method_used = ifelse(method_label == "GMM", ifelse(!is.null(res$method), res$method, "gmm"), res$method))
}

perform_findmarkers_for_celltype <- function(cell_type, method_label) {
  cat("  Analyzing cell type:", cell_type, "with", method_label, "classification...\n")

  cells_to_keep <- which(core_data[[all_cell_types_level]] == cell_type)
  if (length(cells_to_keep) == 0) {
    cat("    No cells found for", cell_type, "\n")
    return(NULL)
  }

  cell_subset <- core_data[, cells_to_keep]
  if (ncol(cell_subset) < min_cells_per_type) {
    cat("    Skipping", cell_type, "- insufficient cells (", ncol(cell_subset), " < ", min_cells_per_type, ")\n", sep = "")
    return(NULL)
  }

  target_expression <- GetAssayData(cell_subset, assay = "RNA", layer = "data")[feature_name, ]
  cls <- classify_cells(target_expression, method_label, cell_type)
  high_expr_cells <- cls$high; low_expr_cells <- cls$low

  if (length(high_expr_cells) < min_cells_per_group || length(low_expr_cells) < min_cells_per_group) {
    cat("    Skipping", cell_type, "- insufficient cells in groups\n")
    return(NULL)
  }

  cell_subset$target_expression_group <- ifelse(colnames(cell_subset) %in% high_expr_cells, "High", "Low")
  Idents(cell_subset) <- "target_expression_group"

  diff_markers <- tryCatch({
    FindMarkers(cell_subset,
      ident.1 = "High", ident.2 = "Low",
      test.use = "wilcox",
      logfc.threshold = logfc_threshold,
      min.pct = min_pct,
      verbose = FALSE)
  }, error = function(e) {
    cat("    Error in FindMarkers:", e$message, "\n"); return(NULL)
  })

  if (is.null(diff_markers)) return(NULL)

  diff_markers$ensembl_id <- rownames(diff_markers)
  if (nrow(gene_mapping) > 0) {
    diff_markers <- diff_markers %>% left_join(gene_mapping, by = "ensembl_id") %>% mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_id, gene_symbol))
  } else {
    diff_markers$gene_symbol <- diff_markers$ensembl_id
  }

  diff_markers <- diff_markers %>% mutate(
    significance = dplyr::case_when(
      p_val_adj <= padj_threshold & avg_log2FC >= logfc_threshold ~ "Upregulated",
      p_val_adj <= padj_threshold & avg_log2FC <= -logfc_threshold ~ "Downregulated",
      TRUE ~ "Not significant"
    ),
    neg_log10_padj = -log10(p_val_adj + 1e-300),
    cell_type = cell_type,
    analysis_method = analysis_method,
    grouping_method = method_label
  )

  return(diff_markers)
}

# =============================================================================
# SECTION 5: RUN ANALYSIS PER METHOD ####
# =============================================================================

run_for_method <- function(method_label) {
  # Per-method output dirs
  out_base <- file.path("step3_differential_analysis")
  ensure_directories(c(
    file.path(out_base, "processed", method_label),
    file.path(out_base, "results", method_label),
    file.path(out_base, "figures", method_label)
  ))

  successful <- c(); failed <- c(); all_diff <- list()

  # Focused first
  for (cell_type in primary_cell_types) {
    if (cell_type %in% all_cell_types) {
      res <- perform_findmarkers_for_celltype(cell_type, method_label)
      if (!is.null(res)) {
        all_diff[[cell_type]] <- res
        successful <- c(successful, cell_type)
        readr::write_csv(res, file.path(out_base, "results", method_label, paste0("findmarkers_", cfg$target_gene_symbol, "_", gsub("[^A-Za-z0-9]", "_", cell_type), "_markers.csv")))
      } else {
        failed <- c(failed, cell_type)
      }
    } else {
      cat("  Cell type", cell_type, "not found in data\n")
      failed <- c(failed, cell_type)
    }
  }

  # Remaining types
  remaining <- setdiff(all_cell_types, primary_cell_types)
  for (cell_type in remaining) {
    res <- perform_findmarkers_for_celltype(cell_type, method_label)
    if (!is.null(res)) {
      all_diff[[cell_type]] <- res
      successful <- c(successful, cell_type)
      readr::write_csv(res, file.path(out_base, "results", method_label, paste0("findmarkers_", cfg$target_gene_symbol, "_", gsub("[^A-Za-z0-9]", "_", cell_type), "_markers.csv")))
    } else {
      failed <- c(failed, cell_type)
    }
  }

  if (length(all_diff) > 0) {
    combined <- dplyr::bind_rows(all_diff)
    readr::write_csv(combined, file.path(out_base, "results", method_label, paste0("findmarkers_", cfg$target_gene_symbol, "_combined_markers.csv")))
  }

  # Summary
  summary_list <- lapply(all_diff, function(df) {
    data.frame(
      cell_type = unique(df$cell_type),
      analysis_method = analysis_method,
      grouping_method = method_label,
      total_genes_tested = nrow(df),
      upregulated_genes = sum(df$significance == "Upregulated"),
      downregulated_genes = sum(df$significance == "Downregulated"),
      not_significant_genes = sum(df$significance == "Not significant"),
      stringsAsFactors = FALSE
    )
  })
  if (length(summary_list) > 0) {
    summary_df <- dplyr::bind_rows(summary_list)
    readr::write_csv(summary_df, file.path(out_base, "results", method_label, paste0("findmarkers_", cfg$target_gene_symbol, "_analysis_summary.csv")))
  }

  invisible(TRUE)
}

for (m in method_order) run_for_method(m)

cat("FindMarkers analysis complete for methods:", paste(method_order, collapse = ", "), "\n")

# =============================================================================
# SECTION 6: JJVOLCANO-STYLE MULTI-PANEL PLOTS ####
# =============================================================================

create_jjvolcano_plot <- function(combined_data, title_suffix = "FindMarkers") {
  # Filter out target gene and prepare data
  plot_data <- combined_data %>%
    dplyr::filter(gene_symbol != cfg$target_gene_symbol & ensembl_id != feature_name) %>%
    dplyr::filter(abs(avg_log2FC) >= logfc_threshold & p_val < 0.05) %>%
    dplyr::mutate(
      type = ifelse(avg_log2FC >= logfc_threshold, "sigUp", "sigDown"),
      cluster = factor(cell_type, levels = primary_cell_types)
    ) %>%
    dplyr::filter(!is.na(cluster))

  if (nrow(plot_data) == 0) return(NULL)

  # Background data range per cluster
  back_data <- purrr::map_dfr(unique(plot_data$cluster), function(x) {
    tmp <- plot_data %>% dplyr::filter(cluster == x)
    if (nrow(tmp) == 0) return(NULL)
    data.frame(cluster = x, min = min(tmp$avg_log2FC) - 0.2, max = max(tmp$avg_log2FC) + 0.2)
  })

  # Top genes per cluster
  top_genes <- plot_data %>% dplyr::group_by(cluster) %>% dplyr::arrange(dplyr::desc(abs(avg_log2FC))) %>% dplyr::slice_head(n = 5) %>% dplyr::ungroup()

  cluster_levels <- levels(plot_data$cluster)
  tile_colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                   "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")[1:length(cluster_levels)]
  names(tile_colors) <- cluster_levels

  p <- ggplot2::ggplot(plot_data, aes(x = cluster, y = avg_log2FC)) +
    ggplot2::geom_col(data = back_data, aes(x = cluster, y = min), fill = "grey93") +
    ggplot2::geom_col(data = back_data, aes(x = cluster, y = max), fill = "grey93") +
    ggplot2::geom_tile(aes(x = cluster, y = 0, fill = cluster), color = "black", height = logfc_threshold * 2, alpha = 0.6, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = tile_colors) +
    ggplot2::geom_jitter(aes(color = type), size = 0.75, width = 0.25, alpha = 0.8) +
    ggplot2::scale_color_manual(values = c("sigDown" = "#0099CC", "sigUp" = "#CC3333"), labels = c("sigDown" = "Down-regulated", "sigUp" = "Up-regulated")) +
    ggrepel::geom_text_repel(data = top_genes, aes(x = cluster, y = avg_log2FC, label = gene_symbol), size = 2.5, max.overlaps = 50, box.padding = 0.3, point.padding = 0.2, force = 3, min.segment.length = 0.05) +
    ggplot2::geom_text(aes(x = cluster, y = 0, label = cluster), size = 3, fontface = "bold", color = "white") +
    ggplot2::scale_y_continuous(n.breaks = 6, expand = expansion(mult = c(0.1, 0.1))) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid = ggplot2::element_blank(),
      legend.position = c(0.85, 0.95), legend.title = ggplot2::element_blank(),
      legend.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.5),
      axis.line.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.5),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16)
    ) +
    ggplot2::labs(title = paste0(title_suffix, ": ", cfg$target_gene_symbol, "-high vs ", cfg$target_gene_symbol, "-low Multi-Panel Comparison"), x = "Clusters", y = "Average log2FoldChange", color = "") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4)))

  return(p)
}

# Generate multi-panel plots per method with high resolution
for (m in method_order) {
  comb_path <- file.path("step3_differential_analysis", "results", m, paste0("findmarkers_", cfg$target_gene_symbol, "_combined_markers.csv"))
  if (file.exists(comb_path)) {
    combined_data <- readr::read_csv(comb_path, show_col_types = FALSE)
    p <- create_jjvolcano_plot(combined_data, "FindMarkers")
    if (!is.null(p)) {
      fig_dir <- file.path("step3_differential_analysis", "figures", m)
      ensure_directories(c(fig_dir))
      pdf_file <- file.path(fig_dir, paste0("findmarkers_jjvolcano_", cfg$target_gene_symbol, "_multipanel.pdf"))
      png_file <- file.path(fig_dir, paste0("findmarkers_jjvolcano_", cfg$target_gene_symbol, "_multipanel.png"))
      ggsave(pdf_file, p, width = 16, height = 10, units = "in")
      ggsave(png_file, p, width = 16, height = 10, units = "in", dpi = 900)
    }
  }
}

cat("High-resolution jjVolcano plots saved for FindMarkers methods where data available.\n")
