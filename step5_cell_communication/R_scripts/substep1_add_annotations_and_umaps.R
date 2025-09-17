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
  library(readr)
  library(ggplot2)
  library(patchwork)
})

source("config/common_utils.R")
cfg <- read_global_config()
set.seed(cfg$seed)

ensure_directories(c(
  "step5_cell_communication/processed",
  "step5_cell_communication/results",
  "step5_cell_communication/figures/dual_umaps",
  "logs"
))

# Determine method order for Step 5 (independent of other steps)
method_setting <- tolower(cfg$step5_group_methods)
method_order <- switch(method_setting,
  gmm = c("GMM"),
  median = c("median"),
  both = if (toupper(cfg$step5_both_primary_method) == "GMM") c("GMM", "median") else c("median", "GMM"),
  c("GMM")
)

# =============================================================================
# SECTION 2: LOAD STEP1 DATA & APPLY CELLCOMM CROSS-TABLE ####
# =============================================================================

cat("Reading core data from Step1...\n")
core_data <- readRDS("step1_add_annotations_from_crosstable/processed/core_data.rds")

cat("Reading cell communication cross table...\n")
cc_path <- cfg$annotation_cross_table_cellchat_path
if (!file.exists(cc_path)) stop("CellChat cross table not found: ", cc_path)
cc_table <- readr::read_csv(cc_path, show_col_types = FALSE)

# Add new annotation level from cross table (unsplit baseline)
base_info <- get_highest_annotation_level(core_data@meta.data)
from_to <- pick_mapping_columns(cc_table, core_data@meta.data)
from_col <- from_to$from
unsplit_col_name <- "annotation_unsplit"
preferred_to <- intersect(c("annotation_cellchat", "annotation_custom", "annotation_target", "annotation_new"), colnames(cc_table))
if (length(preferred_to) > 0) {
  to_col <- preferred_to[1]
} else {
  to_col <- from_to$to
}

cat("Mapping", from_col, "->", to_col, "into", unsplit_col_name, "\n")
core_data <- add_mapped_annotation(core_data, from_col, to_col, cc_table, unsplit_col_name)
new_col_lvl1 <- unsplit_col_name

# =============================================================================
# SECTION 3: SPLIT TARGET CELL TYPES INTO HIGH/LOW PER METHOD ####
# =============================================================================

feature_name <- identify_target_gene_feature(core_data, cfg)
cat("Using target feature:", feature_name, "for symbol:", cfg$target_gene_symbol, "\n")

k_targets <- intersect(cfg$target_cell_types, unique(as.character(core_data@meta.data[[new_col_lvl1]])))
cat("Target cell types for splitting:", paste(k_targets, collapse = ", "), "\n")

current_level_num <- base_info$num + 1
new_cols_created <- c(new_col_lvl1)

add_split_column_for_ct <- function(seurat_obj, match_cell_type_name, gene_symbol, method_label, base_level_col, label_cell_type_name = match_cell_type_name) {
  cells_mask <- which(as.character(seurat_obj@meta.data[[base_level_col]]) == match_cell_type_name)
  if (length(cells_mask) == 0) return(NULL)
  subset_obj <- seurat_obj[, cells_mask]
  expr <- get_expression_vector(subset_obj, feature_name)
  if (method_label == "GMM") res <- split_by_gmm(expr) else res <- split_by_median(expr)
  labels <- label_high_low(colnames(subset_obj), res$high, res$low, label_cell_type_name, gene_symbol)
  list(labels = labels)
}

# Determine which levels to run based on config in a fixed sequence
seq_tokens <- c("unsplit", cfg$target_cell_types, "allsplit")

# Normalize requested tokens
req_tokens <- cfg$step5_levels_to_run
run_unsplit <- any(tolower(req_tokens) == "unsplit")
run_allsplit <- any(tolower(req_tokens) == "allsplit")

normalize_label <- function(x) tolower(gsub("_", "-", x))
avail_cts <- unique(as.character(core_data@meta.data[[new_col_lvl1]]))
avail_norm_to_actual <- setNames(avail_cts, normalize_label(avail_cts))

# Map config target cell types to actual cross-table labels (exact, no M1/M2 remapping)
config_cts <- cfg$target_cell_types
config_to_actual <- vapply(config_cts, function(ct) {
  n <- normalize_label(ct)
  actual <- unname(avail_norm_to_actual[n])
  if (is.na(actual) || is.null(actual)) return(NA_character_) else return(actual)
}, character(1))

# Build ct_sequence preserving config order and availability
ct_sequence <- unique(config_to_actual[!is.na(config_to_actual) & config_to_actual %in% avail_cts])

# Display name map: actual label -> config display label
display_map <- setNames(config_cts, ifelse(is.na(config_to_actual), config_cts, config_to_actual))

# Build selected tokens list strictly from requests, preserving canonical order
canonical_tokens <- c("unsplit", ct_sequence, "allsplit")
requested_norm <- normalize_label(req_tokens)

requested_actuals <- unique(unname(avail_norm_to_actual[requested_norm]))
requested_actuals <- requested_actuals[!is.na(requested_actuals)]

selected_tokens <- character(0)
if ("unsplit" %in% tolower(req_tokens)) selected_tokens <- c(selected_tokens, "unsplit")
for (act in ct_sequence) {
  if (act %in% requested_actuals) selected_tokens <- c(selected_tokens, act)
}
if ("allsplit" %in% tolower(req_tokens)) selected_tokens <- c(selected_tokens, "allsplit")

# Helper to convert actual CT name (from cross-table) back to display CT name from config
display_from_actual <- function(actual_label) {
  if (!is.null(display_map[[actual_label]])) return(display_map[[actual_label]])
  actual_label
}

new_cols_created <- c()
levels_by_method <- list()

# If unsplit requested, create only once and register it
if (any(tolower(selected_tokens) == "unsplit")) {
  new_cols_created <- c(new_cols_created, new_col_lvl1)
  levels_by_method[[length(levels_by_method) + 1]] <- data.frame(level = new_col_lvl1, method = "none", token = "unsplit", stringsAsFactors = FALSE)
}

# Build remaining tokens (excluding unsplit) per method in configured order
tokens_to_build <- selected_tokens[tolower(selected_tokens) != "unsplit"]

build_level_num <- if (any(tolower(selected_tokens) == "unsplit")) (base_info$num + 1) else base_info$num

for (method_label in method_order) {
  for (tok in tokens_to_build) {
    build_level_num <- build_level_num + 1
    # Sanitize column name tokens to avoid hyphens/spaces issues in plotting backends
    if (tolower(tok) == "allsplit") {
      out_col <- paste0("annotation_allsplit_", method_label)
    } else {
      display_ct <- display_from_actual(tok)
      safe_ct <- gsub("[^A-Za-z0-9]+", "_", display_ct)
      out_col <- paste0("annotation_", safe_ct, "_", method_label)
    }
    core_data@meta.data[[out_col]] <- NA_character_

    if (tolower(tok) == "allsplit") {
      # Start from the original cross-table annotation (unchanged for non-target cells)
      core_data@meta.data[[out_col]] <- as.character(core_data@meta.data[[new_col_lvl1]])
      # Overlay splits for each target CT in order, using display names in labels
      for (ct_actual in ct_sequence) {
        display_ct <- display_from_actual(ct_actual)
        res <- add_split_column_for_ct(core_data, ct_actual, cfg$target_gene_symbol, method_label, new_col_lvl1, label_cell_type_name = display_ct)
        if (!is.null(res)) {
          tcells <- names(res$labels)
          core_data@meta.data[tcells, out_col] <- res$labels[tcells]
        }
      }
      core_data@meta.data[[out_col]] <- factor(core_data@meta.data[[out_col]])
      new_cols_created <- c(new_cols_created, out_col)
      levels_by_method[[length(levels_by_method) + 1]] <- data.frame(level = out_col, method = method_label, token = "allsplit", stringsAsFactors = FALSE)
      next
    }

    # Specific CT: split that mapped CT; other cells keep original CT names
  ct_actual <- tok
  display_ct <- display_from_actual(ct_actual)
    core_data@meta.data[[out_col]] <- as.character(core_data@meta.data[[new_col_lvl1]])
    res <- add_split_column_for_ct(core_data, ct_actual, cfg$target_gene_symbol, method_label, new_col_lvl1, label_cell_type_name = display_ct)
    if (!is.null(res)) {
      tcells <- names(res$labels)
      core_data@meta.data[tcells, out_col] <- res$labels[tcells]
    }
    core_data@meta.data[[out_col]] <- factor(core_data@meta.data[[out_col]])
    new_cols_created <- c(new_cols_created, out_col)
    levels_by_method[[length(levels_by_method) + 1]] <- data.frame(level = out_col, method = method_label, token = display_ct, stringsAsFactors = FALSE)
  }
}

# Save updated object with all new levels and level list
saveRDS(core_data, "step5_cell_communication/processed/core_data_cellchat_ready.rds")
readr::write_lines(new_cols_created, "step5_cell_communication/results/new_levels_created.txt")
if (length(levels_by_method) > 0) {
  levels_by_method_df <- dplyr::bind_rows(levels_by_method)
  readr::write_csv(levels_by_method_df, "step5_cell_communication/results/levels_by_method.csv")
}

# =============================================================================
# SECTION 4: DUAL UMAP PLOTS ####
# =============================================================================

create_dual_umap <- function(seurat_obj, annotation1, annotation2, title_suffix) {
  p1 <- Seurat::DimPlot(seurat_obj, group.by = annotation1, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE) +
    ggplot2::ggtitle(annotation1) + ggplot2::theme(legend.position = "right")
  p2 <- Seurat::DimPlot(seurat_obj, group.by = annotation2, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE) +
    ggplot2::ggtitle(annotation2) + ggplot2::theme(legend.position = "right")
  combined <- p1 / p2 + patchwork::plot_annotation(
    title = paste("Dual UMAP:", title_suffix),
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"))
  )
  combined
}

ann3 <- if ("annotation_level_3" %in% colnames(core_data@meta.data)) "annotation_level_3" else NULL

# Determine plot levels with method awareness
levels_df <- tryCatch(readr::read_csv("step5_cell_communication/results/levels_by_method.csv", show_col_types = FALSE), error = function(e) NULL)
plot_entries <- list()
if (!is.null(levels_df) && all(c("level","method") %in% colnames(levels_df))) {
  # Include unsplit once if present
  if ("annotation_unsplit" %in% colnames(core_data@meta.data)) {
    plot_entries[[length(plot_entries) + 1]] <- data.frame(level = "annotation_unsplit", method = "unsplit", stringsAsFactors = FALSE)
  }
  # Add all method-specific levels
  plot_entries[[length(plot_entries) + 1]] <- levels_df %>% dplyr::select(level, method)
  plot_df <- dplyr::bind_rows(plot_entries)
  plot_df <- plot_df %>% dplyr::filter(level %in% colnames(core_data@meta.data))
} else {
  # Fallback: plot any annotation_* columns
  ann_cols <- grep("^annotation_", colnames(core_data@meta.data), value = TRUE)
  plot_df <- data.frame(level = ann_cols, method = "unsplit", stringsAsFactors = FALSE)
}

base_dir <- "step5_cell_communication/figures/dual_umaps"
if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

for (i in seq_len(nrow(plot_df))) {
  lvl <- plot_df$level[i]
  mth <- tolower(plot_df$method[i])
  subdir <- ifelse(mth %in% c("none","unsplit","legacy"), "unsplit", plot_df$method[i])
  out_dir <- file.path(base_dir, subdir)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (!is.null(ann3)) {
    dual <- create_dual_umap(core_data, ann3, lvl, paste(ann3, "vs", lvl))
    out_pdf <- file.path(out_dir, paste0("dual_umap_", ann3, "_vs_", lvl, ".pdf"))
    out_png <- file.path(out_dir, paste0("dual_umap_", ann3, "_vs_", lvl, ".png"))
    ggsave(out_pdf, dual, width = 16, height = 12)
    ggsave(out_png, dual, width = 16, height = 12, dpi = 300)
  }
}

cat("Substep1 complete. Updated core_data, new levels, and dual UMAPs saved.\n")

# NOTE: This unified Substep 1 script replaces the previous consolidated step5_cell_communication.R
# It performs: (1) Apply cellchat cross-table; (2) Create requested split levels per Step 5 config; (3) Save updated object & level list; (4) Generate dual UMAPs.


