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
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(ggsci)
  library(randomcoloR)
})

source("config/common_utils.R")
cfg <- read_global_config()
set.seed(cfg$seed)

# Output dirs
ensure_directories(c(
  "step2_visualization/figures",
  "step2_visualization/results",
  "logs"
))

# =============================================================================
# SECTION 2: READ UPDATED CORE DATA ####
# =============================================================================
cat("Reading updated core data (Step1 output)...\n")
core_data <- readRDS("step1_add_annotations_from_crosstable/processed/core_data.rds")

cat("Core data dimensions:", paste(dim(core_data), collapse = " x "), "\n")
ann_cols <- grep("annotation_level_", colnames(core_data@meta.data), value = TRUE)
cat("Available annotation columns:", paste(ann_cols, collapse = ", "), "\n")

# Detect key annotation columns
ann_info <- get_highest_annotation_level(core_data@meta.data)
ann3 <- if ("annotation_level_3" %in% ann_cols) "annotation_level_3" else ann_cols[pmax(1, length(ann_cols) - 2)]
ann4 <- if ("annotation_level_4" %in% ann_cols) "annotation_level_4" else ann_cols[pmax(1, length(ann_cols) - 1)]
ann5 <- if ("annotation_level_5" %in% ann_cols) "annotation_level_5" else ann_info$col

# =============================================================================
# SECTION 3: VISUALIZE ANNOTATION LEVELS (UMAP) ####
# =============================================================================

create_umap_plot <- function(seurat_obj, group_by, title) {
  seurat_obj@meta.data[[paste0(group_by, "_char")]] <- as.character(seurat_obj@meta.data[[group_by]])
  group_by_char <- paste0(group_by, "_char")
  n_groups <- length(unique(seurat_obj@meta.data[[group_by_char]]))
  if (n_groups <= 16) {
    colors <- c("#AF4034", "#55967e", '#006a8e', "#6a60a9",
                "#D55E00", "#E39E3E", "#0072B2",
                "#CC79A7", '#f9a11b', "#4F86C6", "#fdc23e",
                "#e3632d", "#78C2C4", "#C73E3A", "#2fa1dd", "#f87669")[1:n_groups]
  } else {
    colors <- randomcoloR::distinctColorPalette(n_groups)
  }
  p <- Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = group_by_char,
                 label = TRUE, label.size = 3.5, cols = colors, pt.size = 0.5, raster = FALSE) +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
          legend.position = "right",
          plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold")) +
    ggplot2::ggtitle(title)
  return(p)
}

p1 <- create_umap_plot(core_data, ann3, paste0("Annotation ", ann3))
p2 <- create_umap_plot(core_data, ann4, paste0("Annotation ", ann4))
p3 <- create_umap_plot(core_data, ann5, paste0("Annotation ", ann5))

# Save plots
for (nm in c("annotation_level_3", "annotation_level_4", "annotation_level_5")) {
  if (nm %in% ann_cols) {
    p <- switch(nm, annotation_level_3 = p1, annotation_level_4 = p2, annotation_level_5 = p3)
    ggsave(paste0("step2_visualization/figures/", nm, "_umap.pdf"), plot = p, width = 10, height = 8, dpi = 300)
    ggsave(paste0("step2_visualization/figures/", nm, "_umap.png"), plot = p, width = 10, height = 8, dpi = 300)
  }
}

combined_h <- p1 + p2 + p3 + patchwork::plot_layout(ncol = 3)
if (!is.null(combined_h)) {
  ggsave("step2_visualization/figures/annotation_levels_combined_umap.pdf", plot = combined_h, width = 18, height = 6, dpi = 300)
  ggsave("step2_visualization/figures/annotation_levels_combined_umap.png", plot = combined_h, width = 18, height = 6, dpi = 300)
}

combined_v <- p1 / p3 + patchwork::plot_layout(nrow = 2)
if (!is.null(combined_v)) {
  ggsave("step2_visualization/figures/annotation_level_3_vs_5_vertical_umap.pdf", plot = combined_v, width = 10, height = 16, dpi = 300)
  ggsave("step2_visualization/figures/annotation_level_3_vs_5_vertical_umap.png", plot = combined_v, width = 10, height = 16, dpi = 300)
}

# =============================================================================
# SECTION 4: CELL COUNT BAR PLOTS ####
# =============================================================================

create_cell_count_plot <- function(seurat_obj, group_by, title) {
  cell_counts <- seurat_obj@meta.data %>%
    dplyr::group_by(.data[[group_by]]) %>%
    dplyr::summarise(cell_count = dplyr::n(), .groups = 'drop') %>%
    dplyr::arrange(dplyr::desc(cell_count)) %>%
    dplyr::mutate(percentage = round(cell_count / sum(cell_count) * 100, 1))

  ggplot2::ggplot(cell_counts, aes(x = reorder(.data[[group_by]], cell_count), y = cell_count)) +
    ggplot2::geom_bar(stat = "identity", fill = "#4A90E2", color = "black", linewidth = 0.3) +
    ggplot2::geom_text(aes(label = paste0(cell_count, "\n(", percentage, "%)")), hjust = -0.1, size = 3, color = "black") +
    ggplot2::coord_flip() +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.border = ggplot2::element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid")
    ) +
    ggplot2::labs(x = "Cell Type", y = "Number of Cells", title = title) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
}

plots_counts <- list()
if ("annotation_level_3" %in% ann_cols) {
  plots_counts[["3"]] <- create_cell_count_plot(core_data, "annotation_level_3", "Cell Counts - Annotation Level 3")
  ggsave("step2_visualization/figures/cell_counts_annotation_level_3.pdf", plots_counts[["3"]], width = 12, height = 8, dpi = 300)
  ggsave("step2_visualization/figures/cell_counts_annotation_level_3.png", plots_counts[["3"]], width = 12, height = 8, dpi = 300)
}
if ("annotation_level_4" %in% ann_cols) {
  plots_counts[["4"]] <- create_cell_count_plot(core_data, "annotation_level_4", "Cell Counts - Annotation Level 4")
  ggsave("step2_visualization/figures/cell_counts_annotation_level_4.pdf", plots_counts[["4"]], width = 14, height = 16, dpi = 300)
  ggsave("step2_visualization/figures/cell_counts_annotation_level_4.png", plots_counts[["4"]], width = 14, height = 16, dpi = 300)
}
if (ann5 %in% ann_cols) {
  plots_counts[["5"]] <- create_cell_count_plot(core_data, ann5, paste0("Cell Counts - ", ann5))
  ggsave("step2_visualization/figures/cell_counts_annotation_level_5.pdf", plots_counts[["5"]], width = 12, height = 10, dpi = 300)
  ggsave("step2_visualization/figures/cell_counts_annotation_level_5.png", plots_counts[["5"]], width = 12, height = 10, dpi = 300)
}

if (length(plots_counts) == 3) {
  cell_counts_combined <- plots_counts[["3"]] + plots_counts[["4"]] + plots_counts[["5"]] + patchwork::plot_layout(ncol = 3)
  ggsave("step2_visualization/figures/cell_counts_all_levels_combined.pdf", cell_counts_combined, width = 24, height = 12, dpi = 300)
  ggsave("step2_visualization/figures/cell_counts_all_levels_combined.png", cell_counts_combined, width = 24, height = 12, dpi = 300)
}

# =============================================================================
# SECTION 5: TARGET GENE VISUALIZATIONS ####
# =============================================================================

# Identify target feature robustly
feature_name <- identify_target_gene_feature(core_data, cfg)
cat("Using target feature:", feature_name, "for symbol:", cfg$target_gene_symbol, "\n")

# FeaturePlot (all cells)
p_target_feature <- FeaturePlot(core_data, features = feature_name, reduction = "umap", pt.size = 0.5, raster = FALSE) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  ggtitle(paste0(cfg$target_gene_symbol, " Expression - All Cell Types"))

ggsave(paste0("step2_visualization/figures/", cfg$target_gene_symbol, "_expression_all_cells.pdf"), p_target_feature, width = 10, height = 8, dpi = 300)
sgg <- ggsave(paste0("step2_visualization/figures/", cfg$target_gene_symbol, "_expression_all_cells.png"), p_target_feature, width = 10, height = 8, dpi = 300)

# Violin by ann5
p_target_violin <- VlnPlot(core_data, features = feature_name, group.by = ann5, pt.size = 0) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  ggtitle(paste0(cfg$target_gene_symbol, " Expression Across Cell Types")) +
  NoLegend()

ggsave(paste0("step2_visualization/figures/", cfg$target_gene_symbol, "_violin_all_cells.pdf"), p_target_violin, width = 14, height = 8, dpi = 300)
sgg <- ggsave(paste0("step2_visualization/figures/", cfg$target_gene_symbol, "_violin_all_cells.png"), p_target_violin, width = 14, height = 8, dpi = 300)

# DotPlot by ann5
p_target_dot <- DotPlot(core_data, features = feature_name, group.by = ann5, cols = "RdYlBu") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  ggtitle(paste0(cfg$target_gene_symbol, " Expression - Target Cell Types")) +
  RotatedAxis()

ggsave(paste0("step2_visualization/figures/", cfg$target_gene_symbol, "_dot_all_cells.pdf"), p_target_dot, width = 10, height = 8, dpi = 300)
sgg <- ggsave(paste0("step2_visualization/figures/", cfg$target_gene_symbol, "_dot_all_cells.png"), p_target_dot, width = 10, height = 8, dpi = 300)

# =============================================================================
# SECTION 5.1: TARGETED CELL TYPES SUBSET EXPRESSION (ADDED) ####
# =============================================================================

# Create targeted subset expression visualizations for configured target cell types
available_targets <- intersect(cfg$target_cell_types, unique(as.character(core_data@meta.data[[ann5]])))
if (length(available_targets) > 0) {
  target_mask <- core_data@meta.data[[ann5]] %in% available_targets
  target_cells <- colnames(core_data)[target_mask]
  if (length(target_cells) > 0) {
    core_data_subset <- subset(core_data, cells = target_cells)

    p_target_subset_feature <- FeaturePlot(core_data_subset,
                                           features = feature_name,
                                           reduction = "umap",
                                           pt.size = 0.8,
                                           raster = FALSE) +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5, linetype = "solid"),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
      ggtitle(paste0(cfg$target_gene_symbol, " Expression - Target Cell Types"))

    ggsave(paste0("step2_visualization/figures/", cfg$target_gene_symbol, "_expression_target_cells.pdf"), p_target_subset_feature, width = 10, height = 8, dpi = 300)
    ggsave(paste0("step2_visualization/figures/", cfg$target_gene_symbol, "_expression_target_cells.png"), p_target_subset_feature, width = 10, height = 8, dpi = 300)

    p_target_subset_violin <- VlnPlot(core_data_subset,
                                      features = feature_name,
                                      group.by = ann5,
                                      pt.size = 0) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
      ggtitle(paste0(cfg$target_gene_symbol, " Expression - Target Cell Types")) +
      NoLegend()

    ggsave(paste0("step2_visualization/figures/", cfg$target_gene_symbol, "_violin_target_cells.pdf"), p_target_subset_violin, width = 12, height = 8, dpi = 300)
    ggsave(paste0("step2_visualization/figures/", cfg$target_gene_symbol, "_violin_target_cells.png"), p_target_subset_violin, width = 12, height = 8, dpi = 300)

    p_target_subset_dot <- DotPlot(core_data_subset,
                                   features = feature_name,
                                   group.by = ann5,
                                   cols = "RdYlBu") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
      ggtitle(paste0(cfg$target_gene_symbol, " Expression - Target Cell Types")) +
      RotatedAxis()

    ggsave(paste0("step2_visualization/figures/", cfg$target_gene_symbol, "_dot_target_cells.pdf"), p_target_subset_dot, width = 10, height = 8, dpi = 300)
    ggsave(paste0("step2_visualization/figures/", cfg$target_gene_symbol, "_dot_target_cells.png"), p_target_subset_dot, width = 10, height = 8, dpi = 300)
  }
}

# =============================================================================
# SECTION 6: EXPRESSION PERCENTAGE ACROSS ALL CELL TYPES ####
# =============================================================================

pct_df <- compute_expression_percentage(core_data, feature_name, ann5, threshold = cfg$expression_threshold)
readr::write_csv(pct_df, "step2_visualization/results/expression_percentage_by_celltype.csv")

p_pct <- ggplot2::ggplot(pct_df, aes(x = reorder(.group, percentage), y = percentage)) +
  ggplot2::geom_col(fill = "#4A90E2", color = "black", linewidth = 0.3) +
  ggplot2::coord_flip() +
  ggplot2::theme_classic() +
  ggplot2::theme(
    panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.5),
    axis.title = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
  ) +
  ggplot2::labs(x = "Cell Type", y = paste0("% ", cfg$target_gene_symbol, " Expressed (> ", cfg$expression_threshold, ")"), title = paste0(cfg$target_gene_symbol, " Expression % by Cell Type"))

ggsave("step2_visualization/figures/expression_percentage_by_celltype.pdf", p_pct, width = 10, height = 8)
sgg <- ggsave("step2_visualization/figures/expression_percentage_by_celltype.png", p_pct, width = 10, height = 8, dpi = 300)

# =============================================================================
# SECTION 7: SUMMARY STATS ####
# =============================================================================

annotation_summary <- core_data@meta.data %>%
  group_by(.data[[ann5]]) %>%
  summarise(cell_count = n(), percentage = round(n() / nrow(core_data@meta.data) * 100, 2), .groups = 'drop') %>%
  arrange(desc(cell_count))

readr::write_csv(annotation_summary, "step2_visualization/results/annotation_level_5_summary.csv")

cat("Step 2 visualization complete. Figures and results saved to step2_visualization/\n")
