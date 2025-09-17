# =============================================================================
# Step 4: Pathway Network Visualization - aPEAR Network Plots
# =============================================================================
# Mirrors the referenced implementation. Only loading/saving/naming adapted.
# =============================================================================

# Clear workspace and set working directory
rm(list = ls())
if (getwd() != "/data/y1005/GBmap_Pipeline_v1.0") {
  setwd("/data/y1005/GBmap_Pipeline_v1.0")
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(DOSE)
  library(aPEAR)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(readr)
  library(purrr)
})

# Load config (for consistency and reproducibility seed)
source("config/common_utils.R")
cfg <- read_global_config()
set.seed(cfg$seed)

# GO ontology selection (align with reference). Default to ALL if not set by caller.
if (!exists("GO_ONTOLOGY")) GO_ONTOLOGY <- "ALL"  # BP | MF | CC | ALL
GO_DIR_NAME <- if (toupper(GO_ONTOLOGY) == "ALL") "GO_all" else paste0("GO_", toupper(GO_ONTOLOGY))
GO_DISPLAY_NAME <- if (toupper(GO_ONTOLOGY) == "ALL") "GO-all" else paste0("GO-", toupper(GO_ONTOLOGY))

# Output directories
base_net_dir <- "step4_pathway_enrichment/figures/Network_Plots"
dir.create(base_net_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_net_dir, GO_DIR_NAME), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_net_dir, "KEGG"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_net_dir, "Hallmark"), showWarnings = FALSE, recursive = TRUE)

# 1. Load GSEA Results ####
cat("Loading GSEA enrichment results...\n")

read_if_exists <- function(p) { if (file.exists(p)) readr::read_csv(p, show_col_types = FALSE) else NULL }

combined_go_results <- read_if_exists(paste0("step4_pathway_enrichment/results/combined_", GO_DIR_NAME, "_gsea_results.csv"))
combined_kegg_results <- read_if_exists("step4_pathway_enrichment/results/combined_KEGG_gsea_results.csv")
combined_hallmark_results <- read_if_exists("step4_pathway_enrichment/results/combined_Hallmark_gsea_results.csv")

# Respect config scope for cell types
scope <- tryCatch(tolower(cfg$gsea_celltypes_scope), error = function(e) "focused")
if (!scope %in% c("all", "focused")) scope <- "focused"
target_cts <- cfg$target_cell_types

apply_scope <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if (scope == "focused") {
    df <- df %>% dplyr::filter(cell_type %in% target_cts)
  }
  df
}

combined_go_results <- apply_scope(combined_go_results)
combined_kegg_results <- apply_scope(combined_kegg_results)
combined_hallmark_results <- apply_scope(combined_hallmark_results)

# 2. Helper Functions for Network Visualization ####

prepare_apear_data <- function(gsea_data) {
  if (is.null(gsea_data) || nrow(gsea_data) == 0) return(NULL)
  if ("padj" %in% colnames(gsea_data) && !"p.adjust" %in% colnames(gsea_data)) {
    gsea_data <- gsea_data %>% dplyr::rename(p.adjust = padj)
  }
  significant_data <- gsea_data %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::arrange(p.adjust)
  n_available <- nrow(significant_data)
  n_select <- if (n_available >= 50) 50 else if (n_available >= 40) 40 else if (n_available >= 30) 30 else n_available
  if (n_select < 3) return(NULL)
  apear_data <- significant_data %>%
    dplyr::slice_head(n = n_select) %>%
    dplyr::select(Description, NES, p.adjust, core_enrichment, setSize) %>%
    dplyr::rename(pvalue = p.adjust, geneID = core_enrichment, Count = setSize) %>%
    as.data.frame()
  if (nrow(apear_data) < 3) return(NULL)
  apear_data
}

create_pathway_network <- function(gsea_data, cell_type, method, analysis_type, output_dir) {
  cat("  Creating network plot for", cell_type, "-", method, "-", analysis_type, "\n")
  if (is.null(gsea_data) || nrow(gsea_data) == 0) return(NULL)
  apear_data <- prepare_apear_data(gsea_data)
  if (is.null(apear_data) || nrow(apear_data) < 3) return(NULL)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  clean_cell <- gsub("[^A-Za-z0-9]", "_", cell_type)
  output_prefix <- file.path(output_dir, paste0(clean_cell, "_", method, "_", analysis_type, "_network"))
  res <- list()
  try({
    set.seed(42)
    p1 <- aPEAR::enrichmentNetwork(apear_data, colorBy = "NES", nodeSize = "Count", repelLabels = TRUE, drawEllipses = FALSE)
    ggsave(paste0(output_prefix, "_basic.pdf"), p1, width = 12, height = 10)
    ggsave(paste0(output_prefix, "_basic.png"), p1, width = 12, height = 10, dpi = 300)
    res$basic <- p1
  }, silent = TRUE)
  if (nrow(apear_data) >= 5) {
    tr <- try({
      set.seed(42)
      p2 <- aPEAR::enrichmentNetwork(apear_data, colorBy = "NES", nodeSize = "Count", repelLabels = TRUE, drawEllipses = TRUE)
      ggsave(paste0(output_prefix, "_circled.pdf"), p2, width = 12, height = 10)
      ggsave(paste0(output_prefix, "_circled.png"), p2, width = 12, height = 10, dpi = 300)
      res$circled <- p2
    }, silent = TRUE)
  }
  res
}

create_enhanced_network <- function(gsea_data, cell_type, method, analysis_type, output_dir) {
  apear_data <- prepare_apear_data(gsea_data)
  if (is.null(apear_data) || nrow(apear_data) < 5) return(NULL)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  clean_cell <- gsub("[^A-Za-z0-9]", "_", cell_type)
  output_prefix <- file.path(output_dir, paste0(clean_cell, "_", method, "_", analysis_type, "_enhanced"))
  tryCatch({
    set.seed(654824)
    p3 <- aPEAR::enrichmentNetwork(
      apear_data,
      colorBy = "NES",
      nodeSize = "Count",
      repelLabels = TRUE,
      drawEllipses = TRUE
    ) +
      ggplot2::scale_color_gradientn(
        colours = c("#2166AC", "#4393C3", "#92C5DE", "#F7F7F7", "#FDBF6F", "#E31A1C", "#B30000"),
        name = "NES"
      ) +
      ggplot2::labs(
        title = paste("Enhanced Pathway Network:", cell_type, "-", method, "-", analysis_type),
        subtitle = paste("Number of pathways:", nrow(apear_data))
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5),
        legend.position = "right"
      )
    pdf(paste0(output_prefix, ".pdf"), width = 14, height = 10)
    print(p3)
    dev.off()
    png(paste0(output_prefix, ".png"), width = 14, height = 10, units = "in", res = 300)
    print(p3)
    dev.off()
    p3
  }, error = function(e) NULL)
}

# 3. Generate Network Plots ####

gen_for <- function(df, tag, out_dir) {
  if (is.null(df) || nrow(df) == 0) return(list(n = 0, enh = 0))
  combos <- df %>% dplyr::select(cell_type, method) %>% dplyr::distinct()
  n_ok <- 0
  n_enh <- 0
  for (i in seq_len(nrow(combos))) {
    ct <- combos$cell_type[i]
    m <- combos$method[i]
    cell_method_data <- df %>% dplyr::filter(cell_type == !!ct & method == !!m)
    plots <- create_pathway_network(cell_method_data, ct, m, tag, out_dir)
    if (!is.null(plots)) n_ok <- n_ok + 1
  }
  # Enhanced plots for top combinations by pathway count (cap at 3 as in reference)
  top_combos <- df %>%
    dplyr::group_by(cell_type, method) %>%
    dplyr::summarise(pathway_count = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(pathway_count >= 10) %>%
    dplyr::arrange(dplyr::desc(pathway_count)) %>%
    dplyr::slice_head(n = 3)
  if (nrow(top_combos) > 0) {
    for (i in seq_len(nrow(top_combos))) {
      ct <- top_combos$cell_type[i]
      m <- top_combos$method[i]
      cell_method_data <- df %>% dplyr::filter(cell_type == !!ct & method == !!m)
      p3 <- create_enhanced_network(cell_method_data, ct, m, tag, out_dir)
      if (!is.null(p3)) n_enh <- n_enh + 1
    }
  }
  list(n = n_ok, enh = n_enh)
}

cat("\n=== Generating", GO_DISPLAY_NAME, "Network Plots ===\n")
res_go <- gen_for(combined_go_results, GO_DIR_NAME, file.path(base_net_dir, GO_DIR_NAME))

cat("\n=== Generating KEGG Network Plots ===\n")
res_kegg <- gen_for(combined_kegg_results, "KEGG", file.path(base_net_dir, "KEGG"))

cat("\n=== Generating Hallmark Network Plots ===\n")
res_hall <- gen_for(combined_hallmark_results, "Hallmark", file.path(base_net_dir, "Hallmark"))

# 4. Summary ####

summary_data <- data.frame(
  Analysis_Type = c(GO_DIR_NAME, "KEGG", "Hallmark"),
  Successful_Plot_Sets = c(res_go$n, res_kegg$n, res_hall$n),
  Enhanced_Plots = c(res_go$enh, res_kegg$enh, res_hall$enh),
  Output_Directory = c(
    file.path(base_net_dir, GO_DIR_NAME),
    file.path(base_net_dir, "KEGG"),
    file.path(base_net_dir, "Hallmark")
  )
)
readr::write_csv(summary_data, "step4_pathway_enrichment/results/network_plots_summary.csv")

session_info <- sessionInfo()
saveRDS(session_info, "step4_pathway_enrichment/results/network_plots_session_info.rds")

cat("\n=== Network Plot Generation Complete ===\n")
print(summary_data)
cat("\nOutput directories:\n")
cat("- ", file.path(base_net_dir, GO_DIR_NAME), "\n", sep = "")
cat("- ", file.path(base_net_dir, "KEGG"), "\n", sep = "")
cat("- ", file.path(base_net_dir, "Hallmark"), "\n", sep = "")


