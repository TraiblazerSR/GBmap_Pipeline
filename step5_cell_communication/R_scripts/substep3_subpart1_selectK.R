# Clear workspace and set working directory
rm(list = ls())
if (getwd() != "/data/y1005/GBmap_Pipeline_v1.0") {
  setwd("/data/y1005/GBmap_Pipeline_v1.0")
}

# Initialize logging (match reference style)
LOG_FILE <- "step5_cell_communication/results/selectK_pattern_analysis_log.txt"
log_message <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  entry <- paste0("[", timestamp, "] [", level, "] ", message)
  cat(entry, "\n")
  cat(entry, "\n", file = LOG_FILE, append = TRUE)
}

suppressPackageStartupMessages({
  library(CellChat)
  library(patchwork)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(future)
  if (!require(ggalluvial, quietly = TRUE)) {
    try(install.packages("ggalluvial"), silent = TRUE)
    suppressMessages(library(ggalluvial))
  }
})

source("config/common_utils.R")
cfg <- read_global_config()
set.seed(cfg$seed)

options(device = "pdf")
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 50 * 1024^3)
future::plan("multisession", workers = cfg$workers)

base_dir <- "step5_cell_communication/results"
levels_df <- tryCatch(readr::read_csv(file.path(base_dir, "levels_by_method.csv"), show_col_types = FALSE), error = function(e) NULL)
if (!is.null(levels_df) && all(c("level","method") %in% colnames(levels_df))) {
  iter <- levels_df
} else {
  lv <- readr::read_lines(file.path(base_dir, "new_levels_created.txt"))
  lv <- lv[lv != ""]
  iter <- data.frame(level = lv, method = "legacy", stringsAsFactors = FALSE)
}
if (nrow(iter) == 0) stop("No annotation levels found for selectK")

for (i in seq_len(nrow(iter))) {
  lvl <- iter$level[i]
  mth <- iter$method[i]
  method_dir <- ifelse(tolower(mth) %in% c("none","unsplit","legacy"), "unsplit", mth)
  data_dir <- file.path(base_dir, method_dir, lvl, "data")
  sel_dir <- file.path(base_dir, method_dir, lvl, "selectK_plots")
  if (!dir.exists(sel_dir)) dir.create(sel_dir, recursive = TRUE)
  obj_path <- file.path(data_dir, paste0(lvl, "_cellchat.rds"))
  if (!file.exists(obj_path)) next
  cc <- readRDS(obj_path)
  pathways <- cc@netP$pathways
  if (length(pathways) < 3) next
  try({
    while (dev.cur() > 1) dev.off()
    p_out <- CellChat::selectK(cc, pattern = "outgoing")
    ggsave(file.path(sel_dir, paste0(lvl, "_selectK_outgoing.pdf")), p_out, width = 10, height = 6)
    ggsave(file.path(sel_dir, paste0(lvl, "_selectK_outgoing.png")), p_out, width = 10, height = 6, dpi = 300)
    p_in <- CellChat::selectK(cc, pattern = "incoming")
    ggsave(file.path(sel_dir, paste0(lvl, "_selectK_incoming.pdf")), p_in, width = 10, height = 6)
    ggsave(file.path(sel_dir, paste0(lvl, "_selectK_incoming.png")), p_in, width = 10, height = 6, dpi = 300)
    while (dev.cur() > 1) dev.off()
  }, silent = TRUE)
}

cat("SelectK plots generated. Review and set k in substep3_subpart2_patterns.R\n")


