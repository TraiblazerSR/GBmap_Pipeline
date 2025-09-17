# Clear workspace and set working directory
rm(list = ls())
if (getwd() != "/data/y1005/GBmap_Pipeline_v1.0") {
  setwd("/data/y1005/GBmap_Pipeline_v1.0")
}

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(patchwork)
})

source("config/common_utils.R")
cfg <- read_global_config()
set.seed(cfg$seed)

base_dir <- "step5_cell_communication/results"
levels_df <- tryCatch(readr::read_csv(file.path(base_dir, "levels_by_method.csv"), show_col_types = FALSE), error = function(e) NULL)
if (!is.null(levels_df) && all(c("level","method") %in% colnames(levels_df))) {
  iter <- levels_df
} else {
  lv <- readr::read_lines(file.path(base_dir, "new_levels_created.txt"))
  lv <- lv[lv != ""]
  iter <- data.frame(level = lv, method = "legacy", stringsAsFactors = FALSE)
}
if (nrow(iter) == 0) stop("No annotation levels found for patterns")

# EDIT k values after reviewing selectK plots
k_values <- list(
  # Example defaults; adjust after reviewing selectK plots in each level's selectK_plots/
  # annotation_level_6 = list(outgoing = 3, incoming = 3),
  # annotation_level_7 = list(outgoing = 3, incoming = 3),
  # annotation_level_8 = list(outgoing = 3, incoming = 2),
  # annotation_level_9 = list(outgoing = 3, incoming = 2),
  # annotation_level_10 = list(outgoing = 3, incoming = 2),
  # annotation_level_11 = list(outgoing = 3, incoming = 2),
)

for (i in seq_len(nrow(iter))) {
  lvl <- iter$level[i]
  mth <- iter$method[i]
  method_dir <- ifelse(tolower(mth) %in% c("none","unsplit","legacy"), "unsplit", mth)
  data_dir <- file.path(base_dir, method_dir, lvl, "data")
  pat_dir <- file.path(base_dir, method_dir, lvl, "pattern_analysis")
  if (!dir.exists(pat_dir)) dir.create(pat_dir, recursive = TRUE)
  obj_path <- file.path(data_dir, paste0(lvl, "_cellchat.rds"))
  if (!file.exists(obj_path)) next
  cc <- readRDS(obj_path)
  if (length(cc@netP$pathways) < 3) next

  k_out <- if (!is.null(k_values[[lvl]])) k_values[[lvl]]$outgoing else 3
  k_in  <- if (!is.null(k_values[[lvl]])) k_values[[lvl]]$incoming else 3

  # Outgoing
  if (k_out >= 2) {
    try({
      cc2 <- CellChat::identifyCommunicationPatterns(cc, pattern = "outgoing", k = k_out)
      p_river <- CellChat::netAnalysis_river(cc2, pattern = "outgoing")
      ggsave(file.path(pat_dir, paste0(lvl, "_outgoing_patterns_river.pdf")), p_river, width = 12, height = 14)
      ggsave(file.path(pat_dir, paste0(lvl, "_outgoing_patterns_river.png")), p_river, width = 12, height = 14, dpi = 300)
      p_dot <- CellChat::netAnalysis_dot(cc2, pattern = "outgoing")
      ggsave(file.path(pat_dir, paste0(lvl, "_outgoing_patterns_dot.pdf")), p_dot, width = 18, height = 8)
      ggsave(file.path(pat_dir, paste0(lvl, "_outgoing_patterns_dot.png")), p_dot, width = 18, height = 8, dpi = 300)
      saveRDS(cc2, file.path(data_dir, paste0(lvl, "_cellchat_with_outgoing_patterns.rds")))
    }, silent = TRUE)
  }

  # Incoming
  if (k_in >= 2) {
    try({
      cc3 <- CellChat::identifyCommunicationPatterns(cc, pattern = "incoming", k = k_in)
      p_river <- CellChat::netAnalysis_river(cc3, pattern = "incoming")
      ggsave(file.path(pat_dir, paste0(lvl, "_incoming_patterns_river.pdf")), p_river, width = 12, height = 14)
      ggsave(file.path(pat_dir, paste0(lvl, "_incoming_patterns_river.png")), p_river, width = 12, height = 14, dpi = 300)
      p_dot <- CellChat::netAnalysis_dot(cc3, pattern = "incoming")
      ggsave(file.path(pat_dir, paste0(lvl, "_incoming_patterns_dot.pdf")), p_dot, width = 18, height = 8)
      ggsave(file.path(pat_dir, paste0(lvl, "_incoming_patterns_dot.png")), p_dot, width = 18, height = 8, dpi = 300)
      saveRDS(cc3, file.path(data_dir, paste0(lvl, "_cellchat_with_incoming_patterns.rds")))
    }, silent = TRUE)
  }
}

cat("Pattern analysis completed. Adjust k_values in this script as needed.\n")


