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
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(Matrix)
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
celltype_level <- "annotation_level_5"
all_cell_types <- unique(core_data@meta.data[[celltype_level]])

# Identify target feature robustly
feature_name <- identify_target_gene_feature(core_data, cfg)
cat("Using target feature:", feature_name, "for symbol:", cfg$target_gene_symbol, "\n")

# Gene mapping
all_genes <- rownames(core_data)
tryCatch({
  gene_mapping <- clusterProfiler::bitr(all_genes, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  colnames(gene_mapping) <- c("ensembl_id", "gene_symbol")
}, error = function(e) { gene_mapping <<- data.frame(ensembl_id = character(0), gene_symbol = character(0)) })

# =============================================================================
# SECTION 3: PARAMETERS ####
# =============================================================================

analysis_method <- "Pseudobulk"
min_cells_per_type <- 50
min_cells_per_group <- 20
logfc_threshold <- 0.5
min_samples_per_group <- 3
cells_per_pseudosample <- 50
padj_threshold <- 0.05

primary_cell_types <- cfg$target_cell_types

# =============================================================================
# SECTION 4: CLASSIFICATION & PSEUDOBULK HELPERS ####
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

create_pseudobulk_samples <- function(cell_subset, target_expression, cell_type_name, method_label) {
  cls <- classify_cells(target_expression, method_label, cell_type_name)
  high_expr_cells <- cls$high; low_expr_cells <- cls$low

  if (length(high_expr_cells) < min_cells_per_group || length(low_expr_cells) < min_cells_per_group) return(NULL)

  create_group <- function(cell_names, group_name) {
    if (length(cell_names) < cells_per_pseudosample) return(NULL)
    n_samples <- floor(length(cell_names) / cells_per_pseudosample)
    if (n_samples < min_samples_per_group) return(NULL)
    set.seed(cfg$seed)
    shuffled <- sample(cell_names)
    pb_data <- list(); sample_names <- c()
    for (i in 1:n_samples) {
      s_cells <- shuffled[((i - 1) * cells_per_pseudosample + 1):(i * cells_per_pseudosample)]
      counts <- Matrix::rowSums(GetAssayData(cell_subset, assay = "RNA", layer = "counts")[, s_cells])
      sname <- paste0(group_name, "_sample_", i)
      pb_data[[sname]] <- counts
      sample_names <- c(sample_names, sname)
    }
    list(data = pb_data, samples = sample_names)
  }

  high_pb <- create_group(high_expr_cells, "High")
  low_pb <- create_group(low_expr_cells, "Low")
  if (is.null(high_pb) || is.null(low_pb)) return(NULL)

  all_data <- c(high_pb$data, low_pb$data)
  all_names <- c(high_pb$samples, low_pb$samples)
  mat <- do.call(cbind, all_data)
  colnames(mat) <- all_names

  metadata <- data.frame(sample_id = all_names, group = factor(c(rep("High", length(high_pb$samples)), rep("Low", length(low_pb$samples)))), stringsAsFactors = FALSE)
  rownames(metadata) <- all_names

  list(counts = mat, metadata = metadata, high_samples = length(high_pb$samples), low_samples = length(low_pb$samples), classification_method = cls$method_used)
}

perform_deseq2 <- function(pb) {
  dds <- DESeqDataSetFromMatrix(countData = pb$counts, colData = pb$metadata, design = ~ group)
  min_samples <- min(3, ncol(dds))
  keep <- rowSums(counts(dds) >= 10) >= min_samples
  dds <- dds[keep,]
  dds$group <- relevel(dds$group, ref = "Low")
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group", "High", "Low"))
  as.data.frame(res)
}

perform_pseudobulk_for_celltype <- function(cell_type, method_label) {
  cells_to_keep <- which(core_data[[celltype_level]] == cell_type)
  if (length(cells_to_keep) == 0) return(NULL)
  cell_subset <- core_data[, cells_to_keep]
  if (ncol(cell_subset) < min_cells_per_type) return(NULL)

  target_expression <- GetAssayData(cell_subset, assay = "RNA", layer = "data")[feature_name, ]
  pb <- create_pseudobulk_samples(cell_subset, target_expression, cell_type, method_label)
  if (is.null(pb)) return(NULL)

  res_df <- perform_deseq2(pb)
  if (nrow(res_df) == 0) return(NULL)

  res_df$ensembl_id <- rownames(res_df)
  if (nrow(gene_mapping) > 0) {
    res_df <- res_df %>% dplyr::left_join(gene_mapping, by = "ensembl_id") %>% dplyr::mutate(gene_symbol = ifelse(is.na(gene_symbol), ensembl_id, gene_symbol))
  } else {
    res_df$gene_symbol <- res_df$ensembl_id
  }
  colnames(res_df)[colnames(res_df) == "log2FoldChange"] <- "avg_log2FC"
  colnames(res_df)[colnames(res_df) == "pvalue"] <- "p_val"
  colnames(res_df)[colnames(res_df) == "padj"] <- "p_val_adj"

  res_df$significance <- ifelse(
    res_df$p_val_adj <= padj_threshold & res_df$avg_log2FC >= logfc_threshold, "Upregulated",
    ifelse(res_df$p_val_adj <= padj_threshold & res_df$avg_log2FC <= -logfc_threshold, "Downregulated", "Not significant")
  )
  res_df$neg_log10_padj <- -log10(res_df$p_val_adj + 1e-300)
  res_df$cell_type <- cell_type
  res_df$analysis_method <- analysis_method
  res_df$grouping_method <- method_label

  res_df
}

# =============================================================================
# SECTION 5: RUN ANALYSIS PER METHOD ####
# =============================================================================

run_for_method <- function(method_label) {
  out_base <- file.path("step3_differential_analysis")
  ensure_directories(c(
    file.path(out_base, "processed", method_label),
    file.path(out_base, "results", method_label),
    file.path(out_base, "figures", method_label)
  ))

  all_diff <- list()

  for (cell_type in unique(c(primary_cell_types, setdiff(all_cell_types, primary_cell_types)))) {
    res <- perform_pseudobulk_for_celltype(cell_type, method_label)
    if (!is.null(res)) {
      all_diff[[cell_type]] <- res
      readr::write_csv(res, file.path(out_base, "results", method_label, paste0("pseudobulk_", cfg$target_gene_symbol, "_", gsub("[^A-Za-z0-9]", "_", cell_type), "_markers.csv")))
    }
  }

  if (length(all_diff) > 0) {
    combined <- dplyr::bind_rows(all_diff)
    readr::write_csv(combined, file.path(out_base, "results", method_label, paste0("pseudobulk_", cfg$target_gene_symbol, "_combined_markers.csv")))
  }

  invisible(TRUE)
}

for (m in method_order) run_for_method(m)

cat("Pseudobulk analysis complete for methods:", paste(method_order, collapse = ", "), "\n")

# =============================================================================
# SECTION 6: JJVOLCANO-STYLE MULTI-PANEL PLOTS ####
# =============================================================================

create_jjvolcano_plot <- function(combined_data, title_suffix = "Pseudobulk") {
  plot_data <- combined_data %>%
    dplyr::filter(gene_symbol != cfg$target_gene_symbol & ensembl_id != feature_name) %>%
    dplyr::filter(abs(avg_log2FC) >= logfc_threshold & p_val < 0.05) %>%
    dplyr::mutate(
      type = ifelse(avg_log2FC >= logfc_threshold, "sigUp", "sigDown"),
      cluster = factor(cell_type, levels = primary_cell_types)
    ) %>%
    dplyr::filter(!is.na(cluster))

  if (nrow(plot_data) == 0) return(NULL)

  back_data <- purrr::map_dfr(unique(plot_data$cluster), function(x) {
    tmp <- plot_data %>% dplyr::filter(cluster == x)
    if (nrow(tmp) == 0) return(NULL)
    data.frame(cluster = x, min = min(tmp$avg_log2FC) - 0.2, max = max(tmp$avg_log2FC) + 0.2)
  })

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

for (m in method_order) {
  comb_path <- file.path("step3_differential_analysis", "results", m, paste0("pseudobulk_", cfg$target_gene_symbol, "_combined_markers.csv"))
  if (file.exists(comb_path)) {
    combined_data <- readr::read_csv(comb_path, show_col_types = FALSE)
    p <- create_jjvolcano_plot(combined_data, "Pseudobulk")
    if (!is.null(p)) {
      fig_dir <- file.path("step3_differential_analysis", "figures", m)
      ensure_directories(c(fig_dir))
      pdf_file <- file.path(fig_dir, paste0("pseudobulk_jjvolcano_", cfg$target_gene_symbol, "_multipanel.pdf"))
      png_file <- file.path(fig_dir, paste0("pseudobulk_jjvolcano_", cfg$target_gene_symbol, "_multipanel.png"))
      ggsave(pdf_file, p, width = 16, height = 10, units = "in")
      ggsave(png_file, p, width = 16, height = 10, units = "in", dpi = 900)
    }
  }
}

cat("High-resolution jjVolcano plots saved for Pseudobulk methods where data available.\n")
