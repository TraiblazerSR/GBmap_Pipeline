# Clear workspace and set working directory
rm(list = ls())
if (getwd() != "/data/y1005/GBmap_Pipeline_v1.0") {
  setwd("/data/y1005/GBmap_Pipeline_v1.0")
}

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

source("config/common_utils.R")
cfg <- read_global_config()
set.seed(cfg$seed)

options(device = "pdf")
base_dir <- "step5_cell_communication/results"
levels_df <- tryCatch(readr::read_csv(file.path(base_dir, "levels_by_method.csv"), show_col_types = FALSE), error = function(e) NULL)
if (!is.null(levels_df) && all(c("level","method") %in% colnames(levels_df))) {
  iter <- levels_df
} else {
  lv <- readr::read_lines(file.path(base_dir, "new_levels_created.txt"))
  lv <- lv[lv != ""]
  iter <- data.frame(level = lv, method = "legacy", stringsAsFactors = FALSE)
}
if (nrow(iter) == 0) stop("No annotation levels found for bubble plots")

safe_device_close <- function() {
  tryCatch({ while (dev.cur() > 1) dev.off() }, error = function(e) {})
}

# Move stray root-level png/svg into the level's bubble_plots dir
cleanup_stray_files <- function(output_dir, level_name) {
  tryCatch({
    stray <- list.files(getwd(), pattern = "\\.(png|svg)$", full.names = TRUE)
    if (length(stray) > 0) {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      for (f in stray) {
        file.rename(f, file.path(output_dir, paste0(level_name, "_stray_", basename(f))))
      }
    }
  }, error = function(e) {})
}

safe_bubble_plot <- function(cellchat_obj, output_path_base, plot_params) {
  tryCatch({
    safe_device_close()
    p <- CellChat::netVisual_bubble(
      cellchat_obj,
      sources.use = plot_params$sources.use,
      targets.use = plot_params$targets.use,
      signaling = plot_params$signaling,
      pairLR.use = plot_params$pairLR.use,
      remove.isolate = plot_params$remove.isolate,
      sort.by.source = plot_params$sort.by.source,
      sort.by.target = plot_params$sort.by.target,
      font.size = plot_params$font.size,
      font.size.title = plot_params$font.size.title
    )
    ggplot2::ggsave(paste0(output_path_base, ".pdf"), p, width = plot_params$width, height = plot_params$height)
    ggplot2::ggsave(paste0(output_path_base, ".png"), p, width = plot_params$width, height = plot_params$height, dpi = 300)
    TRUE
  }, error = function(e) { FALSE })
}

for (i in seq_len(nrow(iter))) {
  lvl <- iter$level[i]
  mth <- iter$method[i]
  method_dir <- ifelse(tolower(mth) %in% c("none","unsplit","legacy"), "unsplit", mth)
  data_dir <- file.path(base_dir, method_dir, lvl, "data")
  bubble_dir <- file.path(base_dir, method_dir, lvl, "bubble_plots")
  if (!dir.exists(bubble_dir)) dir.create(bubble_dir, recursive = TRUE)
  cleanup_stray_files(bubble_dir, lvl)
  obj_path <- file.path(data_dir, paste0(lvl, "_cellchat.rds"))
  if (!file.exists(obj_path)) next
  cc <- readRDS(obj_path)
  groups <- levels(cc@idents)
  if (length(groups) < 2) next
  pathways <- cc@netP$pathways
  src_idx <- 1:min(2, length(groups))
  tgt_idx <- max(1, min(3, length(groups))):length(groups)

  params <- list(
    sources.use = src_idx, targets.use = tgt_idx, signaling = NULL, pairLR.use = NULL,
    remove.isolate = FALSE, sort.by.source = FALSE, sort.by.target = FALSE,
    font.size = 10, font.size.title = 12, width = 14, height = 15
  )
  safe_bubble_plot(cc, file.path(bubble_dir, paste0(lvl, "_bubble_interactions")), params)

  if (length(pathways) >= 2) {
    params$signaling <- pathways[1:min(3, length(pathways))]
    safe_bubble_plot(cc, file.path(bubble_dir, paste0(lvl, "_bubble_pathways")), params)
  }

  # L-R pairs
  pairLR <- NULL
  if (length(pathways) >= 1) {
    pairLR <- tryCatch(CellChat::extractEnrichedLR(cc, signaling = pathways[1:min(2, length(pathways))], geneLR.return = FALSE), error = function(e) NULL)
    if (!is.null(pairLR) && nrow(pairLR) > 20) pairLR <- pairLR[1:20, , drop = FALSE]
  }
  if (!is.null(pairLR) && nrow(pairLR) > 0) {
    params$signaling <- NULL
    params$pairLR.use <- pairLR
    params$remove.isolate <- TRUE
    params$font.size <- 9
    params$font.size.title <- 11
    safe_bubble_plot(cc, file.path(bubble_dir, paste0(lvl, "_bubble_LR_pairs")), params)

    params$sources.use <- NULL
    params$targets.use <- NULL
    params$sort.by.target <- TRUE
    safe_bubble_plot(cc, file.path(bubble_dir, paste0(lvl, "_bubble_sort_by_target")), params)

    params$sources.use <- src_idx
    params$sort.by.target <- FALSE
    params$sort.by.source <- TRUE
    safe_bubble_plot(cc, file.path(bubble_dir, paste0(lvl, "_bubble_sort_by_source")), params)

    if (length(groups) >= 4) {
      params$targets.use <- tgt_idx
      params$sort.by.target <- TRUE
      params$font.size <- 8
      params$font.size.title <- 10
      safe_bubble_plot(cc, file.path(bubble_dir, paste0(lvl, "_bubble_sort_by_both")), params)
    }
  }

  if (length(pathways) >= 5) {
    params <- list(
      sources.use = NULL, targets.use = NULL, signaling = pathways[1:min(8, length(pathways))], pairLR.use = NULL,
      remove.isolate = FALSE, sort.by.source = FALSE, sort.by.target = FALSE,
      font.size = 8, font.size.title = 10, width = 16, height = 18
    )
    safe_bubble_plot(cc, file.path(bubble_dir, paste0(lvl, "_bubble_pathways_overview")), params)
  }
}

cat("Bubble plots generated for all levels.\n")


