# Clear workspace and set working directory
rm(list = ls())
if (getwd() != "/data/y1005/GBmap_Pipeline_v1.0") {
  setwd("/data/y1005/GBmap_Pipeline_v1.0")
}

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

source("config/common_utils.R")
cfg <- read_global_config()
set.seed(cfg$seed)

options(device = "pdf")
options(stringsAsFactors = FALSE)

# Move any stray root-level png/svg plots into this level's plots dir
cleanup_stray_files <- function(output_dir, level_name) {
  tryCatch({
    stray <- list.files(getwd(), pattern = "\\.(png|svg)$", full.names = TRUE)
    if (length(stray) > 0) {
      plots_dir <- file.path(output_dir, "plots")
      if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
      for (f in stray) {
        file.rename(f, file.path(plots_dir, paste0(level_name, "_stray_", basename(f))))
      }
    }
  }, error = function(e) {})
}

# Levels to analyze (method-aware)
levels_by_method_path <- "step5_cell_communication/results/levels_by_method.csv"
levels_df <- tryCatch(readr::read_csv(levels_by_method_path, show_col_types = FALSE), error = function(e) NULL)

if (is.null(levels_df) || !all(c("level", "method") %in% colnames(levels_df))) {
  legacy_levels <- readr::read_lines("step5_cell_communication/results/new_levels_created.txt")
  legacy_levels <- legacy_levels[nzchar(legacy_levels)]
  stopifnot(length(legacy_levels) > 0)
  levels_df <- data.frame(level = legacy_levels, method = "legacy", stringsAsFactors = FALSE)
}

plot_level_lr_pairs <- function(level_name, method_label) {
  method_dir <- ifelse(tolower(method_label) %in% c("none", "unsplit", "legacy"), "unsplit", method_label)
  base_results_dir <- file.path("step5_cell_communication/results", method_dir)
  data_dir <- file.path(base_results_dir, level_name, "data")
  lr_dir <- file.path(base_results_dir, level_name, "LR_pair_analysis")
  ensure_directories(c(lr_dir))

  cc_path <- file.path(data_dir, paste0(level_name, "_cellchat.rds"))
  if (!file.exists(cc_path)) {
    message("Skip ", level_name, " [", method_label, "]: missing ", cc_path)
    return(invisible(FALSE))
  }
  cc <- readRDS(cc_path)

  pathways <- cc@netP$pathways
  if (length(pathways) == 0) return(invisible(TRUE))

  vertex.receiver <- seq(1, min(4, length(levels(cc@idents))))
  plots_dir <- file.path(base_results_dir, level_name, "plots")
  ensure_directories(c(plots_dir))
  cleanup_stray_files(file.path(base_results_dir, level_name), level_name)
  # Batch visualizations for all pathways (use pdf for non-gg; ggsave only for gg-returning)
  batch_dir <- file.path(base_results_dir, level_name, "batch_pathway_analysis")
  ensure_directories(c(batch_dir))
  for (pathway_name in pathways) {
    # Hierarchy (non-gg)
    pdf(file.path(batch_dir, paste0(level_name, "_", pathway_name, "_batch_hierarchy.pdf")), width = 10, height = 8)
    try(CellChat::netVisual(cc, signaling = pathway_name, vertex.receiver = vertex.receiver, layout = "hierarchy"), silent = TRUE)
    dev.off()
    # Circle (non-gg)
    pdf(file.path(batch_dir, paste0(level_name, "_", pathway_name, "_batch_circle.pdf")), width = 8, height = 8)
    try(CellChat::netVisual(cc, signaling = pathway_name, vertex.receiver = vertex.receiver, layout = "circle"), silent = TRUE)
    dev.off()
    # LR contribution (gg)
    p_contrib <- try(CellChat::netAnalysis_contribution(cc, signaling = pathway_name), silent = TRUE)
    if (!inherits(p_contrib, "try-error")) {
      ggsave(file.path(batch_dir, paste0(level_name, "_", pathway_name, "_batch_LR_contribution.pdf")), p_contrib, width = 8, height = 6)
      ggsave(file.path(batch_dir, paste0(level_name, "_", pathway_name, "_batch_LR_contribution.png")), p_contrib, width = 8, height = 6, dpi = 300)
    }
  }
  for (pathway in pathways) {
    pairLR <- tryCatch(CellChat::extractEnrichedLR(cc, signaling = pathway, geneLR.return = FALSE), error = function(e) NULL)
    if (is.null(pairLR) || nrow(pairLR) == 0) next
    for (i in seq_len(nrow(pairLR))) {
      LR.show <- pairLR[i, ]
      lr_tag <- paste0(LR.show$ligand, "_", LR.show$receptor)
      pdf(file.path(lr_dir, paste0(level_name, "_", pathway, "_", lr_tag, "_individual_hierarchy.pdf")), width = 10, height = 8)
      try(CellChat::netVisual_individual(cc, signaling = pathway, pairLR.use = LR.show, vertex.receiver = vertex.receiver, layout = "hierarchy"), silent = TRUE)
      dev.off()
      pdf(file.path(lr_dir, paste0(level_name, "_", pathway, "_", lr_tag, "_individual_circle.pdf")), width = 8, height = 8)
      try(CellChat::netVisual_individual(cc, signaling = pathway, pairLR.use = LR.show, layout = "circle"), silent = TRUE)
      dev.off()
      pdf(file.path(lr_dir, paste0(level_name, "_", pathway, "_", lr_tag, "_individual_chord.pdf")), width = 8, height = 8)
      try(CellChat::netVisual_individual(cc, signaling = pathway, pairLR.use = LR.show, layout = "chord"), silent = TRUE)
      dev.off()
    }
  }
  # Chord gene analysis (non-gg)
  try({
    n_groups <- length(levels(cc@idents))
    if (n_groups >= 4) {
      sources_use <- 1:min(2, n_groups - 2)
      targets_use <- (min(3, n_groups - 1)):n_groups
      pdf(file.path(plots_dir, paste0(level_name, "_chord_gene_interactions.pdf")), width = 10, height = 8)
      CellChat::netVisual_chord_gene(cc, sources.use = sources_use, targets.use = targets_use, lab.cex = 0.5)
      dev.off()
    }
  }, silent = TRUE)
  # Multi-group chord (non-gg)
  try({
    cell_groups <- levels(cc@idents)
    group_mapping <- rep("Group1", length(cell_groups))
    if (length(cell_groups) >= 6) {
      group_mapping[1:(length(cell_groups)%/%3)] <- "Group1"
      group_mapping[((length(cell_groups)%/%3)+1):(2*(length(cell_groups)%/%3))] <- "Group2"
      group_mapping[((2*(length(cell_groups)%/%3))+1):length(cell_groups)] <- "Group3"
    }
    names(group_mapping) <- cell_groups
    if (length(pathways) > 0) {
      pdf(file.path(plots_dir, paste0(level_name, "_multigroup_chord.pdf")), width = 10, height = 8)
      CellChat::netVisual_chord_cell(cc, signaling = pathways[1], group = group_mapping, title.name = paste0(pathways[1], " signaling network"))
      dev.off()
    }
  }, silent = TRUE)
  cleanup_stray_files(file.path(base_results_dir, level_name), level_name)
  invisible(TRUE)
}

for (i in seq_len(nrow(levels_df))) {
  lvl <- levels_df$level[i]
  mth <- levels_df$method[i]
  try(plot_level_lr_pairs(lvl, mth), silent = TRUE)
}

cat("Substep 5 (batch + individual LR pairs) complete. All pathways and pairs plotted.\n")


