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
  library(Matrix)
  library(readr)
})

source("config/common_utils.R")
cfg <- read_global_config()
set.seed(cfg$seed)

# Output directories
ensure_directories(c(
  "step0_data_exploration/results",
  "step0_data_exploration/figures",
  "logs"
))

# =============================================================================
# SECTION 2: LOAD CORE DATA ####
# =============================================================================
cat("=== Loading core.rds data ===\n")
core_data <- readRDS(cfg$core_data_path)

# Basic info
cat("Class:", class(core_data)[1], "\n")
cat("Cells:", ncol(core_data), " Genes:", nrow(core_data), "\n")

# =============================================================================
# SECTION 3: STRUCTURE SUMMARY ####
# =============================================================================

sink("step0_data_exploration/results/core_data_structure_summary.txt")
cat("=== SEURAT OBJECT STRUCTURE SUMMARY ===\n")
cat("Generated on:", as.character(Sys.time()), "\n\n")
cat("Class:", class(core_data)[1], "\n")
cat("Size:", format(object.size(core_data), units = "GB"), "\n")
cat("Cells:", ncol(core_data), "\n")
cat("Features:", nrow(core_data), "\n")
cat("Assays:", paste(names(core_data@assays), collapse = ", "), "\n")
cat("Reductions:", paste(names(core_data@reductions), collapse = ", "), "\n")
cat("Meta.data columns:", ncol(core_data@meta.data), "\n")
sink()

# =============================================================================
# SECTION 4: DETAILED STRUCTURE (ADAPTED) ####
# =============================================================================

create_detailed_structure <- function(obj, obj_name = "core_data", indent = "", max_depth = 10, current_depth = 0) {
  if (current_depth > max_depth) return(character(0))
  result <- character(0)
  if (current_depth == 0) {
    obj_class <- class(obj)[1]
    obj_size <- format(object.size(obj), units = "GB")
    result <- c(result, paste0(obj_name, " Large ", obj_class, " (", obj_size, ")"))
  }
  if (isS4(obj)) {
    slots <- slotNames(obj)
    for (slot_name in slots) {
      slot_obj <- slot(obj, slot_name)
      if (slot_name == "assays") {
        result <- c(result, paste0(indent, "..@assays :List of ", length(slot_obj)))
        assay_names <- names(slot_obj)
        for (assay_name in assay_names) {
          assay_obj <- slot_obj[[assay_name]]
          assay_class <- class(assay_obj)[1]
          result <- c(result, paste0(indent, ".. ..$ ", assay_name, ": Formal class '", assay_class, "' with ", length(slotNames(assay_obj)), " slots"))
        }
      } else if (slot_name == "meta.data") {
        result <- c(result, paste0(indent, "..@meta.data :'data.frame': ", nrow(slot_obj), " obs. of ", ncol(slot_obj), " variables:"))
        for (i in 1:min(25, ncol(slot_obj))) {
          var_name <- names(slot_obj)[i]
          var_class <- class(slot_obj[[var_name]])[1]
          result <- c(result, paste0(indent, ".. ..$ ", var_name, " : ", var_class))
        }
      } else if (slot_name == "reductions") {
        result <- c(result, paste0(indent, "..@reductions :List of ", length(slot_obj)))
        if (length(slot_obj) > 0) {
          for (rn in names(slot_obj)) {
            red_obj <- slot_obj[[rn]]
            result <- c(result, paste0(indent, ".. ..$ ", rn, ": Formal class '", class(red_obj)[1], "' with ", length(slotNames(red_obj)), " slots"))
          }
        }
      } else {
        # Print compact line for other slots
        result <- c(result, paste0(indent, "..@", slot_name, ": ", class(slot_obj)[1]))
      }
    }
  }
  result
}

detailed_structure <- create_detailed_structure(core_data)
readr::write_lines(detailed_structure, "step0_data_exploration/results/detailed_core_structure.txt")

cat("Structure summaries saved under step0_data_exploration/results\n")
