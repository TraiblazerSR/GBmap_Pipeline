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
})

source("config/common_utils.R")
cfg <- read_global_config()
set.seed(cfg$seed)

ensure_directories(c(
  "step1_add_annotations_from_crosstable/processed",
  "step1_add_annotations_from_crosstable/results",
  "step1_add_annotations_from_crosstable/figures",
  "logs"
))

# =============================================================================
# SECTION 2: LOAD DATA & CROSS-TABLE ####
# =============================================================================
cat("Reading core data...\n")
core_data <- readRDS(cfg$core_data_path)

cat("Reading annotation cross table for Step1...\n")
cross_path <- cfg$annotation_cross_table_step1_path
if (!file.exists(cross_path)) stop("Cross table not found: ", cross_path)
annotation_cross_table <- readr::read_csv(cross_path, show_col_types = FALSE)

# =============================================================================
# SECTION 3: CREATE NEXT ANNOTATION LEVEL VIA MAPPING ####
# =============================================================================
info <- get_highest_annotation_level(core_data@meta.data)
cat("Highest existing annotation level:", ifelse(is.null(info$col), "None", paste0(info$col, " (", info$num, ")")), "\n")

mapping_cols <- pick_mapping_columns(annotation_cross_table, core_data@meta.data)
from_col <- mapping_cols$from
# The new level must be the next numeric level
new_col <- paste0("annotation_level_", info$num + 1)

# Prefer these 'to' columns if present, otherwise fall back to detection
preferred_to <- intersect(c("annotation_custom", "annotation_cellchat", "annotation_TAMs", "annotation_target", "annotation_new"), colnames(annotation_cross_table))
if (length(preferred_to) > 0) {
  to_col <- preferred_to[1]
} else {
  to_col <- mapping_cols$to
}

cat("Mapping from:", from_col, " -> ", to_col, " into ", new_col, "\n")

core_data <- add_mapped_annotation(core_data, from_col, to_col, annotation_cross_table, new_col)

# Verification summary
cat("Annotation added:", new_col, "\n")
summary_tbl <- core_data@meta.data %>% dplyr::group_by(.data[[new_col]]) %>% dplyr::summarise(cells = dplyr::n(), .groups = 'drop') %>% dplyr::arrange(dplyr::desc(cells))
readr::write_csv(summary_tbl, file.path("step1_add_annotations_from_crosstable/results", paste0(new_col, "_summary.csv")))

# Save updated object
saveRDS(core_data, "step1_add_annotations_from_crosstable/processed/core_data.rds")

cat("Saved updated core_data with ", new_col, " to step1_add_annotations_from_crosstable/processed/core_data.rds\n")
