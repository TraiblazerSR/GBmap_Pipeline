# Clear workspace and set working directory
rm(list = ls())
if (getwd() != "/data/y1005/GBmap_Pipeline_v1.0") {
  setwd("/data/y1005/GBmap_Pipeline_v1.0")
}

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(future)
  if (!require(ggalluvial, quietly = TRUE)) {
    cat("Installing ggalluvial package...\n")
    try(install.packages("ggalluvial"), silent = TRUE)
    suppressWarnings(suppressMessages(library(ggalluvial)))
  }
  if (!require(reticulate, quietly = TRUE)) {
    cat("Note: reticulate not available - UMAP similarity analysis may be skipped\n")
  }
})

source("config/common_utils.R")
cfg <- read_global_config()
set.seed(cfg$seed)

# Logging
LOG_FILE <- "step5_cell_communication/results/cellchat_analysis_log.txt"
log_message <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  entry <- paste0("[", timestamp, "] [", level, "] ", message)
  cat(entry, "\n")
  cat(entry, "\n", file = LOG_FILE, append = TRUE)
}

# Device behavior
options(device = "pdf")
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 50 * 1024^3)
future::plan("multisession", workers = cfg$workers)

# Ensure base dirs
ensure_directories(c(
  "step5_cell_communication/results",
  "step5_cell_communication/results/logs"
))

# Load Seurat object
cat("Loading updated Seurat object with new annotation levels...\n")
seurat_path <- "step5_cell_communication/processed/core_data_cellchat_ready.rds"
if (!file.exists(seurat_path)) stop("Seurat object not found for CellChat: ", seurat_path)
core_data <- readRDS(seurat_path)

# Load DB and PPI
data(CellChatDB.human)
CellChatDB.use <- CellChatDB.human
data(PPI.human)

# Utility: cleanup stray images
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

# Convert features to symbols where possible (matches reference behavior)
convert_to_gene_symbols <- function(seurat_obj) {
  obj <- seurat_obj
  if ("feature_name" %in% colnames(obj@assays$RNA@meta.features)) {
    gene_symbols <- as.character(obj@assays$RNA@meta.features$feature_name)
    names(gene_symbols) <- rownames(obj@assays$RNA@meta.features)
    valid <- !is.na(gene_symbols) & gene_symbols != "" & gene_symbols != "NA"
    gene_symbols <- gene_symbols[valid]
    unique_ok <- !duplicated(gene_symbols)
    gene_symbols <- gene_symbols[unique_ok]
    obj2 <- obj[names(gene_symbols), ]
    obj2@assays$RNA@meta.features$original_id <- rownames(obj2@assays$RNA@meta.features)
    rownames(obj2@assays$RNA@meta.features) <- gene_symbols
    try({
      if (class(obj2@assays$RNA)[1] == "Assay5") {
        dimnames(obj2@assays$RNA@layers$counts)[[1]] <- gene_symbols
        dimnames(obj2@assays$RNA@layers$data)[[1]] <- gene_symbols
      } else {
        dimnames(obj2@assays$RNA@counts)[[1]] <- gene_symbols
        dimnames(obj2@assays$RNA@data)[[1]] <- gene_symbols
      }
    }, silent = TRUE)
    return(obj2)
  }
  obj
}

create_output_dirs <- function(base_dir, level_name) {
  level_dir <- file.path(base_dir, level_name)
  subdirs <- c("data", "plots", "reports", "pathways", "systems_analysis", "LR_pair_analysis", "chord_gene_plots", "batch_pathway_analysis")
  ensure_directories(c(level_dir, file.path(level_dir, subdirs)))
  level_dir
}

# Levels to analyze (prefer levels_by_method.csv for method-aware outputs)
levels_by_method_path <- "step5_cell_communication/results/levels_by_method.csv"
levels_df <- tryCatch(readr::read_csv(levels_by_method_path, show_col_types = FALSE), error = function(e) NULL)

if (!is.null(levels_df) && all(c("level", "method") %in% colnames(levels_df))) {
  levels_df <- levels_df %>% dplyr::filter(level %in% colnames(core_data@meta.data))
  if (nrow(levels_df) == 0) stop("No valid levels found in levels_by_method.csv")
  cat("Annotation levels (method-aware):", paste(paste0(levels_df$level, "[", levels_df$method, "]"), collapse = ", "), "\n")
} else {
  # Fallback to legacy list
  legacy_levels <- readr::read_lines("step5_cell_communication/results/new_levels_created.txt")
  legacy_levels <- legacy_levels[legacy_levels %in% colnames(core_data@meta.data)]
  if (length(legacy_levels) == 0) stop("No new annotation levels recorded to analyze.")
  levels_df <- data.frame(level = legacy_levels, method = "legacy", stringsAsFactors = FALSE)
  cat("Annotation levels (legacy):", paste(levels_df$level, collapse = ", "), "\n")
}

run_level_analysis <- function(level_name, method_label) {
  method_dir <- ifelse(tolower(method_label) %in% c("none", "unsplit", "legacy"), "unsplit", method_label)
  base_results_dir <- file.path("step5_cell_communication/results", method_dir)
  output_dir <- create_output_dirs(base_results_dir, level_name)
  # Convert to symbols for CellChat compatibility
  core_data_for_cellchat <- convert_to_gene_symbols(core_data)
  # Create CellChat object
  cc <- tryCatch({
    CellChat::createCellChat(object = core_data_for_cellchat, group.by = level_name, assay = "RNA")
  }, error = function(e) {
    log_message(paste("Error creating CellChat object for", level_name, ":", e$message), "ERROR")
    return(NULL)
  })
  if (is.null(cc)) return(NULL)
  cc@DB <- CellChatDB.use
  # Preprocess
  cc <- CellChat::subsetData(cc)
  cc <- tryCatch(CellChat::identifyOverExpressedGenes(cc), error = function(e) cc)
  cc <- tryCatch(CellChat::identifyOverExpressedInteractions(cc), error = function(e) cc)
  cc <- tryCatch(CellChat::smoothData(cc, adj = PPI.human), error = function(e) cc)
  cc <- tryCatch(CellChat::computeCommunProb(cc, raw.use = TRUE), error = function(e) cc)
  cc <- tryCatch(CellChat::filterCommunication(cc, min.cells = 5), error = function(e) cc)
  cc <- tryCatch(CellChat::computeCommunProbPathway(cc), error = function(e) cc)
  cc <- tryCatch(CellChat::aggregateNet(cc), error = function(e) cc)
  # Save object and CSV
  data_dir <- file.path(output_dir, "data")
  saveRDS(cc, file.path(data_dir, paste0(level_name, "_cellchat.rds")))
  comm <- tryCatch(CellChat::subsetCommunication(cc), error = function(e) NULL)
  if (!is.null(comm) && nrow(comm) > 0) {
    write.csv(comm, file.path(data_dir, paste0(level_name, "_communications.csv")), row.names = FALSE)
  }
  # Basic network visualizations (non-gg; keep pdf)
  while (dev.cur() > 1) dev.off()
  groupSize <- as.numeric(table(cc@idents))
  pdf(file.path(output_dir, "plots", paste0(level_name, "_network_circle_plots.pdf")), width = 16, height = 8)
  par(mfrow = c(1,2), xpd = TRUE)
  CellChat::netVisual_circle(cc@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
  CellChat::netVisual_circle(cc@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")
  dev.off()
  # Individual signaling patterns
  while (dev.cur() > 1) dev.off()
  pdf(file.path(output_dir, "plots", paste0(level_name, "_individual_signaling.pdf")), width = 20, height = 15)
  mat <- cc@net$weight
  n_groups <- nrow(mat)
  par(mfrow = c(ceiling(n_groups/4), 4), xpd = TRUE)
  for (i in seq_len(n_groups)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    CellChat::netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  dev.off()
  # Pathway-specific visualizations
  pathways <- cc@netP$pathways
  pathway_dir <- file.path(output_dir, "pathways")
  for (pathway in pathways) {
    try({
      while (dev.cur() > 1) dev.off()
      vertex.receiver <- seq(1, min(4, length(levels(cc@idents))))
      pdf(file.path(pathway_dir, paste0(level_name, "_", pathway, "_hierarchy.pdf")), width = 10, height = 8)
      CellChat::netVisual_aggregate(cc, signaling = pathway, vertex.receiver = vertex.receiver, layout = "hierarchy")
      dev.off()
      pdf(file.path(pathway_dir, paste0(level_name, "_", pathway, "_circle.pdf")), width = 8, height = 8)
      CellChat::netVisual_aggregate(cc, signaling = pathway, layout = "circle")
      dev.off()
      try({
        pdf(file.path(pathway_dir, paste0(level_name, "_", pathway, "_chord.pdf")), width = 10, height = 10)
        CellChat::netVisual_aggregate(cc, signaling = pathway, layout = "chord")
        dev.off()
      }, silent = TRUE)
      pdf(file.path(pathway_dir, paste0(level_name, "_", pathway, "_heatmap.pdf")), width = 8, height = 6)
      CellChat::netVisual_heatmap(cc, signaling = pathway, color.heatmap = "Reds")
      dev.off()
      p_contrib <- CellChat::netAnalysis_contribution(cc, signaling = pathway)
      ggsave(file.path(pathway_dir, paste0(level_name, "_", pathway, "_LR_contribution.pdf")), p_contrib, width = 8, height = 6)
      ggsave(file.path(pathway_dir, paste0(level_name, "_", pathway, "_LR_contribution.png")), p_contrib, width = 8, height = 6, dpi = 300)
      while (dev.cur() > 1) dev.off()
    }, silent = TRUE)
  }
  # Chord gene analysis and Multi-group chord moved to Substep5 (substep5_batch_and_individual_lr_pairs.R)
  # Individual L-R pair visualizations moved to Substep5 (substep5_batch_and_individual_lr_pairs.R)
  # Gene lists instead of expression plots
  if (length(pathways) > 0) {
    for (pathway in pathways[1:min(3, length(pathways))]) {
      try({
        genes_file <- file.path(output_dir, "data", paste0(level_name, "_", pathway, "_genes.txt"))
        genes <- CellChat::extractEnrichedLR(cc, signaling = pathway, geneLR.return = TRUE)
        if (length(genes) > 0) writeLines(genes, genes_file)
      }, silent = TRUE)
    }
  }
  # Batch processing moved to Substep5 (substep5_batch_and_individual_lr_pairs.R)
  # Systems analysis
  systems_dir <- file.path(output_dir, "systems_analysis")
  try({
    cc <- CellChat::netAnalysis_computeCentrality(cc, slot.name = "netP")
    if (length(pathways) > 0) {
      pdf(file.path(systems_dir, paste0(level_name, "_centrality_heatmap.pdf")), width = 12, height = 6)
      CellChat::netAnalysis_signalingRole_network(cc, signaling = pathways[1:min(5, length(pathways))], width = 12, height = 6, font.size = 10)
      dev.off()
    }
    gg1 <- CellChat::netAnalysis_signalingRole_scatter(cc)
    if (length(pathways) >= 2) {
      gg2 <- CellChat::netAnalysis_signalingRole_scatter(cc, signaling = pathways[1:2])
      p_scatter <- gg1 + gg2
    } else {
      p_scatter <- gg1
    }
    ggsave(file.path(systems_dir, paste0(level_name, "_signaling_roles_scatter.pdf")), p_scatter, width = 16, height = 8)
    ggsave(file.path(systems_dir, paste0(level_name, "_signaling_roles_scatter.png")), p_scatter, width = 16, height = 8, dpi = 300)
    pdf(file.path(systems_dir, paste0(level_name, "_signaling_roles_heatmap.pdf")), width = 20, height = 28)
    par(cex = 0.7, cex.axis = 0.7, cex.lab = 0.8, cex.main = 0.9)
    ht1 <- CellChat::netAnalysis_signalingRole_heatmap(cc, pattern = "outgoing", font.size = 10, font.size.title = 12, width = 10, height = 28)
    ht2 <- CellChat::netAnalysis_signalingRole_heatmap(cc, pattern = "incoming", font.size = 10, font.size.title = 12, width = 10, height = 28)
    print(ht1 + ht2)
    dev.off()
  }, silent = TRUE)
  # Similarity analysis with UMAP (commented out per request)
  # try({
  #   umap_available <- tryCatch({ reticulate::py_module_available("umap") }, error = function(e) FALSE)
  #   if (umap_available) {
  #     while (dev.cur() > 1) dev.off()
  #     cc <- CellChat::computeNetSimilarity(cc, type = "functional")
  #     cc <- CellChat::netEmbedding(cc, type = "functional")
  #     cc <- CellChat::netClustering(cc, type = "functional")
  #     p_fun <- CellChat::netVisual_embedding(cc, type = "functional", label.size = 3.5)
  #     ggsave(file.path(systems_dir, paste0(level_name, "_functional_similarity.pdf")), p_fun, width = 8, height = 6)
  #     ggsave(file.path(systems_dir, paste0(level_name, "_functional_similarity.png")), p_fun, width = 8, height = 6, dpi = 300)
  #     cc <- CellChat::computeNetSimilarity(cc, type = "structural")
  #     cc <- CellChat::netEmbedding(cc, type = "structural")
  #     cc <- CellChat::netClustering(cc, type = "structural")
  #     p_str <- CellChat::netVisual_embedding(cc, type = "structural", label.size = 3.5)
  #     ggsave(file.path(systems_dir, paste0(level_name, "_structural_similarity.pdf")), p_str, width = 8, height = 6)
  #     ggsave(file.path(systems_dir, paste0(level_name, "_structural_similarity.png")), p_str, width = 8, height = 6, dpi = 300)
  #     p_fun_zoom <- CellChat::netVisual_embeddingZoomIn(cc, type = "functional", nCol = 2)
  #     ggsave(file.path(systems_dir, paste0(level_name, "_functional_similarity_zoomin.pdf")), p_fun_zoom, width = 12, height = 8)
  #     ggsave(file.path(systems_dir, paste0(level_name, "_functional_similarity_zoomin.png")), p_fun_zoom, width = 12, height = 8, dpi = 300)
  #     p_str_zoom <- CellChat::netVisual_embeddingZoomIn(cc, type = "structural", nCol = 2)
  #     ggsave(file.path(systems_dir, paste0(level_name, "_structural_similarity_zoomin.pdf")), p_str_zoom, width = 12, height = 8)
  #     ggsave(file.path(systems_dir, paste0(level_name, "_structural_similarity_zoomin.png")), p_str_zoom, width = 12, height = 8, dpi = 300)
  #   }
  # }, silent = TRUE)
  # Save final object and cleanup
  saveRDS(cc, file.path(data_dir, paste0(level_name, "_cellchat.rds")))
  cleanup_stray_files(output_dir, level_name)
  while (dev.cur() > 1) dev.off()
  TRUE
}

for (i in seq_len(nrow(levels_df))) {
  lvl <- levels_df$level[i]
  mth <- levels_df$method[i]
  log_message(paste("Starting CellChat analysis for", lvl, "[", mth, "]"))
  ok <- run_level_analysis(lvl, mth)
  if (isTRUE(ok)) log_message(paste("Completed CellChat analysis for", lvl, "[", mth, "]")) else log_message(paste("Skipped", lvl, "[", mth, "]"), "WARN")
}

cat("Substep 2 complete. Full CellChat analysis and figures saved.\n")


