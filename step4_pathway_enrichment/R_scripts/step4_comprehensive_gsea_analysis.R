# ================================================================================
# Step 4: Comprehensive GSEA Enrichment Analysis - Focused Cell Types
# ================================================================================
# This script performs GSEA enrichment analysis on DEG results from Step 3
# for focused cell types and mirrors the referenced implementation. Only
# loading, saving, naming, and working directory are adapted to this project.
# ================================================================================
# Clear workspace and set working directory
rm(list = ls())
project_name = "GBmap_Pipeline_v1.0"
if (getwd() != "/data/y1005/GBmap_Pipeline_v1.0") {
  setwd("/data/y1005/GBmap_Pipeline_v1.0")
}

# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(msigdbr)
library(fgsea)
library(ggrepel)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)

# Load global config and seed
source("config/common_utils.R")
cfg <- read_global_config()
set.seed(cfg$seed)

# *** CONFIGURABLE GO ONTOLOGY ***
# Define GO ontology (same behavior as reference; adjust manually here if needed)
if (!exists("GO_ONTOLOGY")) GO_ONTOLOGY <- "ALL"  # Options: "BP", "MF", "CC", "ALL"

# Dynamic GO ontology variables
GO_DIR_NAME <- if(GO_ONTOLOGY == "ALL") "GO_all" else paste0("GO_", GO_ONTOLOGY)
GO_DISPLAY_NAME <- if(GO_ONTOLOGY == "ALL") "GO-all" else paste0("GO-", GO_ONTOLOGY)

# Create output directories
if (!dir.exists("step4_pathway_enrichment/results")) {
  dir.create("step4_pathway_enrichment/results", recursive = TRUE)
}
if (!dir.exists("step4_pathway_enrichment/figures")) {
  dir.create("step4_pathway_enrichment/figures", recursive = TRUE)
}
if (!dir.exists(paste0("step4_pathway_enrichment/figures/", GO_DIR_NAME))) {
  dir.create(paste0("step4_pathway_enrichment/figures/", GO_DIR_NAME), recursive = TRUE)
}
if (!dir.exists("step4_pathway_enrichment/figures/KEGG")) {
  dir.create("step4_pathway_enrichment/figures/KEGG", recursive = TRUE)
}
if (!dir.exists("step4_pathway_enrichment/figures/Focused_Pathways")) {
  dir.create("step4_pathway_enrichment/figures/Focused_Pathways", recursive = TRUE)
}
if (!dir.exists("step4_pathway_enrichment/figures/Hallmark")) {
  dir.create("step4_pathway_enrichment/figures/Hallmark", recursive = TRUE)
}

# ================================================================================
# 1. Load and Prepare DEG Data ####
cat("Loading DEG results from Step 3...\n")

# Helper to load combined results (adapts loading only; preserves reference semantics)
load_combined <- function(prefix) {
  # Prefer legacy single combined file at root if present
    legacy <- file.path("step3_differential_analysis/results", paste0(prefix, "_", cfg$target_gene_symbol, "_combined_markers.csv"))
  if (file.exists(legacy)) return(readr::read_csv(legacy, show_col_types = FALSE))
  # Otherwise, bind across method subdirectories (GMM/median)
  method_dirs <- list.dirs("step3_differential_analysis/results", full.names = TRUE, recursive = FALSE)
    files <- unlist(lapply(method_dirs, function(md) Sys.glob(file.path(md, paste0(prefix, "_", cfg$target_gene_symbol, "_combined_markers.csv")))))
  if (length(files) == 0) stop("No combined markers files found for prefix ", prefix)
  dplyr::bind_rows(lapply(files, readr::read_csv, show_col_types = FALSE))
}

# Load combined DEG results (pseudobulk and findmarkers)
pseudobulk_results <- load_combined("pseudobulk")
findmarkers_results <- load_combined("findmarkers")

# Function to prepare gene list for GSEA
prepare_gene_list <- function(deg_data, cell_type = NULL) {
  if (!is.null(cell_type)) {
    deg_data <- deg_data %>% dplyr::filter(cell_type == !!cell_type)
    cat("    Filtered data for", cell_type, ":", nrow(deg_data), "genes\n")
  }
  gene_list <- deg_data$avg_log2FC
  names(gene_list) <- deg_data$gene_symbol
  gene_list <- gene_list[!duplicated(names(gene_list))]
  gene_list <- sort(gene_list, decreasing = TRUE)
  return(gene_list)
}

# *** FOCUSED CELL TYPES (configurable) ***
focused_cell_types <- if (!is.null(cfg$target_cell_types)) cfg$target_cell_types else c("Endothelial", "TAM-BDM", "TAM-MG")

# Get available cell types from both datasets that match focused types
cell_types_pseudobulk <- unique(pseudobulk_results$cell_type)
cell_types_findmarkers <- unique(findmarkers_results$cell_type)
all_available_types <- unique(c(cell_types_pseudobulk, cell_types_findmarkers))

# Determine scope from config (all|focused)
scope <- tryCatch(tolower(cfg$gsea_celltypes_scope), error = function(e) "focused")
if (!scope %in% c("all", "focused")) scope <- "focused"
all_cell_types <- if (scope == "all") all_available_types else intersect(focused_cell_types, all_available_types)

cat("Focused cell types defined (from config):", paste(focused_cell_types, collapse = ", "), "\n")
cat("Available cell types from data:", paste(all_available_types, collapse = ", "), "\n")
cat("GSEA scope:", scope, "\n")
cat("Final cell types for analysis (", length(all_cell_types), "):", paste(all_cell_types, collapse = ", "), "\n")

# ================================================================================
# 2. Load Focused Pathways ####
cat("Loading focused pathways...\n")

focused_pathways <- readr::read_csv(cfg$focused_pathways_path, col_names = c("Category", "Pathway"))
focused_pathways <- focused_pathways %>% 
  dplyr::filter(!is.na(Category) & !is.na(Pathway) & Category != "" & Pathway != "")

msigdb_hallmarks <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")

convert_pathway_name <- function(pathway_name) {
  pathway_name %>%
    stringr::str_to_upper() %>%
    stringr::str_replace_all("[^A-Z0-9]", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_remove("^_|_$")
}

focused_pathways$MSigDB_Name <- paste0("HALLMARK_", sapply(focused_pathways$Pathway, convert_pathway_name))

focused_msigdb <- msigdb_hallmarks %>%
  dplyr::filter(gs_name %in% focused_pathways$MSigDB_Name)

cat("Found", length(unique(focused_msigdb$gs_name)), "matching MSigDB hallmark pathways\n")

# ================================================================================
# 3. GSEA Analysis Functions ####

# Function to perform GO GSEA
perform_go_gsea <- function(gene_list, cell_type_name, method_name) {
  cat("Performing", GO_DISPLAY_NAME, "GSEA for", cell_type_name, "using", method_name, "...\n")
  tryCatch({
    gsea_result <- clusterProfiler::gseGO(
      geneList = gene_list,
      OrgDb = org.Hs.eg.db,
      ont = GO_ONTOLOGY,
      keyType = "SYMBOL",
      minGSSize = 15,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      verbose = FALSE
    )
    if (nrow(gsea_result@result) > 0) {
      gsea_result@result$cell_type <- cell_type_name
      gsea_result@result$method <- method_name
      return(gsea_result@result)
    } else {
      return(NULL)
    }
  }, error = function(e) {
    cat("Error in", GO_DISPLAY_NAME, "GSEA for", cell_type_name, ":", e$message, "\n")
    return(NULL)
  })
}

# Function to perform KEGG GSEA
perform_kegg_gsea <- function(gene_list, cell_type_name, method_name) {
  cat("Performing KEGG GSEA for", cell_type_name, "using", method_name, "...\n")
  tryCatch({
    gene_df <- data.frame(SYMBOL = names(gene_list), FC = gene_list)
    gene_df$ENTREZID <- clusterProfiler::bitr(gene_df$SYMBOL, 
                                              fromType = "SYMBOL", 
                                              toType = "ENTREZID", 
                                              OrgDb = org.Hs.eg.db)$ENTREZID[match(gene_df$SYMBOL, 
                                                                                   clusterProfiler::bitr(gene_df$SYMBOL, 
                                                                                                         fromType = "SYMBOL", 
                                                                                                         toType = "ENTREZID", 
                                                                                                         OrgDb = org.Hs.eg.db)$SYMBOL)]
    gene_df <- gene_df[!is.na(gene_df$ENTREZID), ]
    entrez_list <- gene_df$FC
    names(entrez_list) <- gene_df$ENTREZID
    entrez_list <- sort(entrez_list, decreasing = TRUE)
    gsea_result <- clusterProfiler::gseKEGG(
      geneList = entrez_list,
      organism = "hsa",
      minGSSize = 15,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      verbose = FALSE
    )
    if (nrow(gsea_result@result) > 0) {
      gsea_result@result$cell_type <- cell_type_name
      gsea_result@result$method <- method_name
      return(gsea_result@result)
    } else {
      return(NULL)
    }
  }, error = function(e) {
    cat("Error in KEGG GSEA for", cell_type_name, ":", e$message, "\n")
    return(NULL)
  })
}

# Function to perform focused MSigDB GSEA
perform_focused_gsea <- function(gene_list, cell_type_name, method_name) {
  cat("Performing Focused MSigDB GSEA for", cell_type_name, "using", method_name, "...\n")
  tryCatch({
    focused_gene_sets <- split(focused_msigdb$gene_symbol, focused_msigdb$gs_name)
    fgsea_result <- fgsea::fgsea(
      pathways = focused_gene_sets,
      stats = gene_list,
      minSize = 15,
      maxSize = 500,
      eps = 0
    )
    if (nrow(fgsea_result) > 0) {
      fgsea_result <- fgsea_result %>%
        dplyr::mutate(
          Description = pathway,
          p.adjust = padj,
          setSize = size,
          core_enrichment = sapply(leadingEdge, paste, collapse = "/")
        )
      fgsea_result$cell_type <- cell_type_name
      fgsea_result$method <- method_name
      fgsea_result$category <- focused_pathways$Category[match(fgsea_result$Description, focused_pathways$MSigDB_Name)]
      fgsea_result <- fgsea_result %>% dplyr::select(-pathway, -padj, -size, -leadingEdge)
      return(as.data.frame(fgsea_result))
    } else {
      return(NULL)
    }
  }, error = function(e) {
    cat("Error in Focused GSEA for", cell_type_name, ":", e$message, "\n")
    return(NULL)
  })
}

# Function to perform Hallmark GSEA
perform_hallmark_gsea <- function(gene_list, cell_type_name, method_name) {
  cat("Performing Hallmark GSEA for", cell_type_name, "using", method_name, "...\n")
  tryCatch({
    hallmark_gene_sets <- split(msigdb_hallmarks$gene_symbol, msigdb_hallmarks$gs_name)
    fgsea_result <- fgsea::fgsea(
      pathways = hallmark_gene_sets,
      stats = gene_list,
      minSize = 15,
      maxSize = 500,
      eps = 0
    )
    if (nrow(fgsea_result) > 0) {
      fgsea_result <- fgsea_result %>%
        dplyr::mutate(
          Description = pathway,
          p.adjust = padj,
          setSize = size,
          core_enrichment = sapply(leadingEdge, paste, collapse = "/")
        ) %>%
        dplyr::select(-pathway, -padj, -size, -leadingEdge)
      fgsea_result$cell_type <- cell_type_name
      fgsea_result$method <- method_name
      return(as.data.frame(fgsea_result))
    } else {
      return(NULL)
    }
  }, error = function(e) {
    cat("Error in Hallmark GSEA for", cell_type_name, ":", e$message, "\n")
    return(NULL)
  })
}

# ================================================================================
# 4. Run GSEA Analysis for All Cell Types and Methods ####

# Initialize result storage
all_go_results <- list()
all_kegg_results <- list()
all_focused_results <- list()
all_hallmark_results <- list()

# Process Pseudobulk results for focused cell types only
cat("\n=== Processing Pseudobulk Results ===\n")
focused_pseudobulk_types <- intersect(all_cell_types, cell_types_pseudobulk)
cat("Processing", length(focused_pseudobulk_types), "focused cell types from pseudobulk data\n")
cat("Focused pseudobulk types:", paste(focused_pseudobulk_types, collapse = ", "), "\n")

for (cell_type in focused_pseudobulk_types) {
  cat("  Processing pseudobulk cell type:", cell_type, "\n")
  gene_list <- prepare_gene_list(pseudobulk_results, cell_type)
  if (length(gene_list) > 0) {
    go_result <- perform_go_gsea(gene_list, cell_type, "Pseudobulk")
    if (!is.null(go_result)) {
      all_go_results[[paste(cell_type, "Pseudobulk", sep = "_")]] <- go_result
    }
    kegg_result <- perform_kegg_gsea(gene_list, cell_type, "Pseudobulk")
    if (!is.null(kegg_result)) {
      all_kegg_results[[paste(cell_type, "Pseudobulk", sep = "_")]] <- kegg_result
    }
    focused_result <- perform_focused_gsea(gene_list, cell_type, "Pseudobulk")
    if (!is.null(focused_result)) {
      all_focused_results[[paste(cell_type, "Pseudobulk", sep = "_")]] <- focused_result
    }
    hallmark_result <- perform_hallmark_gsea(gene_list, cell_type, "Pseudobulk")
    if (!is.null(hallmark_result)) {
      all_hallmark_results[[paste(cell_type, "Pseudobulk", sep = "_")]] <- hallmark_result
    }
  }
}

# Process FindMarkers results for focused cell types only
cat("\n=== Processing FindMarkers Results ===\n")
focused_findmarkers_types <- intersect(all_cell_types, cell_types_findmarkers)
cat("Processing", length(focused_findmarkers_types), "focused cell types from findmarkers data\n")
cat("Focused findmarkers types:", paste(focused_findmarkers_types, collapse = ", "), "\n")

for (cell_type in focused_findmarkers_types) {
  cat("  Processing findmarkers cell type:", cell_type, "\n")
  gene_list <- prepare_gene_list(findmarkers_results, cell_type)
  if (length(gene_list) > 0) {
    go_result <- perform_go_gsea(gene_list, cell_type, "FindMarkers")
    if (!is.null(go_result)) {
      all_go_results[[paste(cell_type, "FindMarkers", sep = "_")]] <- go_result
    }
    kegg_result <- perform_kegg_gsea(gene_list, cell_type, "FindMarkers")
    if (!is.null(kegg_result)) {
      all_kegg_results[[paste(cell_type, "FindMarkers", sep = "_")]] <- kegg_result
    }
    focused_result <- perform_focused_gsea(gene_list, cell_type, "FindMarkers")
    if (!is.null(focused_result)) {
      all_focused_results[[paste(cell_type, "FindMarkers", sep = "_")]] <- focused_result
    }
    hallmark_result <- perform_hallmark_gsea(gene_list, cell_type, "FindMarkers")
    if (!is.null(hallmark_result)) {
      all_hallmark_results[[paste(cell_type, "FindMarkers", sep = "_")]] <- hallmark_result
    }
  }
}

# ================================================================================
# 5. Combine Results ####

combined_go_results <- do.call(rbind, all_go_results)
combined_kegg_results <- do.call(rbind, all_kegg_results)
combined_focused_results <- do.call(rbind, all_focused_results)
combined_hallmark_results <- do.call(rbind, all_hallmark_results)

cat("\nGSEA Analysis Summary:\n")
if (!is.null(combined_go_results)) {
  cat(GO_DISPLAY_NAME, "results:", nrow(combined_go_results), "enriched pathways\n")
} else {
  cat(GO_DISPLAY_NAME, "results: 0 enriched pathways\n")
}
if (!is.null(combined_kegg_results)) {
  cat("KEGG results:", nrow(combined_kegg_results), "enriched pathways\n")
} else {
  cat("KEGG results: 0 enriched pathways\n")
}
if (!is.null(combined_focused_results)) {
  cat("Focused results:", nrow(combined_focused_results), "enriched pathways\n")
} else {
  cat("Focused results: 0 enriched pathways\n")
}
if (!is.null(combined_hallmark_results)) {
  cat("Hallmark results:", nrow(combined_hallmark_results), "enriched pathways\n")
} else {
  cat("Hallmark results: 0 enriched pathways\n")
}

# Save combined results
readr::write_csv(combined_go_results, paste0("step4_pathway_enrichment/results/combined_", GO_DIR_NAME, "_gsea_results.csv"))
readr::write_csv(combined_kegg_results, "step4_pathway_enrichment/results/combined_KEGG_gsea_results.csv")
readr::write_csv(combined_focused_results, "step4_pathway_enrichment/results/combined_focused_gsea_results.csv")
readr::write_csv(combined_hallmark_results, "step4_pathway_enrichment/results/combined_Hallmark_gsea_results.csv")

# ================================================================================
# 6. Visualization Functions ####

create_lollipop_plot_by_nes <- function(gsea_data, title_prefix, output_path) {
  if (nrow(gsea_data) == 0) return(NULL)
  if ("padj" %in% colnames(gsea_data) && !"p.adjust" %in% colnames(gsea_data)) {
    gsea_data <- gsea_data %>% dplyr::rename(p.adjust = padj)
  }
  plot_data <- gsea_data %>%
    dplyr::mutate(abs_NES = abs(NES)) %>%
    dplyr::arrange(desc(abs_NES)) %>%
    dplyr::slice_head(n = 20) %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(
      Description = stringr::str_wrap(Description, width = 50),
      Description = paste0(Description, " (", row_number(), ")"),
      Description = factor(Description, levels = Description)
    )
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = NES, y = Description)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, y = Description, yend = Description), color = "grey60", linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(color = -log10(p.adjust), size = setSize), alpha = 0.8) +
    ggplot2::scale_color_viridis_c(name = "-log10(padj)", option = "plasma") +
    ggplot2::scale_size_continuous(name = "Gene Count", range = c(2, 8)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
    ggplot2::labs(title = paste(title_prefix, "- Top 20 Pathways by |NES|"), x = "Normalized Enrichment Score (NES)", y = "Pathway") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 9),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )
  ggplot2::ggsave(paste0(output_path, "_lollipop_by_NES.pdf"), p, width = 12, height = 8)
  ggplot2::ggsave(paste0(output_path, "_lollipop_by_NES.png"), p, width = 12, height = 8, dpi = 300)
  return(p)
}

create_lollipop_plot_by_significance <- function(gsea_data, title_prefix, output_path) {
  if (nrow(gsea_data) == 0) return(NULL)
  if ("padj" %in% colnames(gsea_data) && !"p.adjust" %in% colnames(gsea_data)) {
    gsea_data <- gsea_data %>% dplyr::rename(p.adjust = padj)
  }
  plot_data <- gsea_data %>%
    dplyr::mutate(neg_log10_padj = -log10(p.adjust)) %>%
    dplyr::arrange(desc(neg_log10_padj)) %>%
    dplyr::slice_head(n = 20) %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(
      Description = stringr::str_wrap(Description, width = 50),
      Description = paste0(Description, " (", row_number(), ")"),
      Description = factor(Description, levels = Description)
    )
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = NES, y = Description)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, y = Description, yend = Description), color = "grey60", linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(color = neg_log10_padj, size = setSize), alpha = 0.8) +
    ggplot2::scale_color_viridis_c(name = "-log10(padj)", option = "plasma") +
    ggplot2::scale_size_continuous(name = "Gene Count", range = c(2, 8)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
    ggplot2::labs(title = paste(title_prefix, "- Top 20 Pathways by Significance"), x = "Normalized Enrichment Score (NES)", y = "Pathway") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 9),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )
  ggplot2::ggsave(paste0(output_path, "_lollipop_by_significance.pdf"), p, width = 12, height = 8)
  ggplot2::ggsave(paste0(output_path, "_lollipop_by_significance.png"), p, width = 12, height = 8, dpi = 300)
  return(p)
}

create_focused_pathway_heatmap <- function(focused_data, title_prefix, output_path) {
  if (nrow(focused_data) == 0) return(NULL)
  if ("p.adjust" %in% colnames(focused_data) && !"padj" %in% colnames(focused_data)) {
    focused_data <- focused_data %>% dplyr::rename(padj = p.adjust)
  }
  focused_pathways_order <- readr::read_csv(cfg$focused_pathways_path, col_names = c("Category", "Pathway"))
  focused_pathways_order <- focused_pathways_order %>% 
    dplyr::filter(!is.na(Category) & !is.na(Pathway) & Category != "" & Pathway != "") %>%
    dplyr::mutate(
      MSigDB_Name = paste0("HALLMARK_", sapply(Pathway, function(x) {
        x %>%
          stringr::str_to_upper() %>%
          stringr::str_replace_all("[^A-Z0-9]", "_") %>%
          stringr::str_replace_all("_+", "_") %>%
          stringr::str_remove("^_|_$")
      })),
      pathway_clean = stringr::str_to_title(stringr::str_replace_all(Pathway, "_", " "))
    )
  plot_data <- focused_data %>%
    dplyr::select(Description, NES, padj, cell_type, method, category) %>%
    dplyr::mutate(
      comparison = paste(cell_type, method, sep = "_"),
      neg_log10_padj = -log10(padj),
      pathway_clean = stringr::str_remove(Description, "HALLMARK_"),
      pathway_clean = stringr::str_replace_all(pathway_clean, "_", " "),
      pathway_clean = stringr::str_to_title(pathway_clean)
    ) %>%
    dplyr::left_join(
      focused_pathways_order %>% dplyr::select(MSigDB_Name, Category) %>% dplyr::rename(Description = MSigDB_Name, category_ordered = Category),
      by = "Description"
    ) %>%
    dplyr::mutate(category = ifelse(!is.na(category_ordered), category_ordered, category)) %>%
    dplyr::select(-category_ordered)
  pathway_order <- focused_pathways_order %>% dplyr::filter(pathway_clean %in% plot_data$pathway_clean) %>% dplyr::pull(pathway_clean)
  additional_pathways <- setdiff(unique(plot_data$pathway_clean), pathway_order)
  pathway_order <- c(pathway_order, additional_pathways)
  nes_matrix <- plot_data %>% dplyr::select(pathway_clean, comparison, NES) %>% tidyr::pivot_wider(names_from = comparison, values_from = NES, values_fill = 0) %>% tibble::column_to_rownames("pathway_clean") %>% as.matrix()
  padj_matrix <- plot_data %>% dplyr::select(pathway_clean, comparison, neg_log10_padj) %>% tidyr::pivot_wider(names_from = comparison, values_from = neg_log10_padj, values_fill = 0) %>% tibble::column_to_rownames("pathway_clean") %>% as.matrix()
  common_pathways <- intersect(rownames(nes_matrix), rownames(padj_matrix))
  common_comparisons <- intersect(colnames(nes_matrix), colnames(padj_matrix))
  available_pathways <- intersect(pathway_order, common_pathways)
  nes_matrix <- nes_matrix[available_pathways, common_comparisons, drop = FALSE]
  padj_matrix <- padj_matrix[available_pathways, common_comparisons, drop = FALSE]
  if (nrow(nes_matrix) == 0 || ncol(nes_matrix) == 0) return(NULL)
  nes_colors <- circlize::colorRamp2(c(min(nes_matrix, na.rm = TRUE), 0, max(nes_matrix, na.rm = TRUE)), c("blue", "grey90", "red"))
  n_categories <- length(unique(plot_data$category))
  category_colors <- if (n_categories <= 8) RColorBrewer::brewer.pal(max(3, n_categories), "Set3")[1:n_categories] else rainbow(n_categories)
  names(category_colors) <- unique(plot_data$category)
  size_fun <- function(x) { pmax(0.2, pmin(1, x / max(padj_matrix, na.rm = TRUE))) }
  size_legend <- ComplexHeatmap::Legend(title = "-log10(padj)", title_gp = grid::gpar(fontsize = 10), labels = as.character(c(1,2,3,5)), labels_gp = grid::gpar(fontsize = 9), type = "points", pch = 16, size = grid::unit(size_fun(c(1,2,3,5)) * 10, "mm"), legend_height = grid::unit(4, "cm"))
  pdf(paste0(output_path, "_focused_heatmap.pdf"), width = 14, height = 10)
  ht <- ComplexHeatmap::Heatmap(
    nes_matrix,
    name = "NES",
    col = nes_colors,
    rect_gp = grid::gpar(type = "none"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.circle(
        x = x, y = y,
        r = size_fun(padj_matrix[i, j]) * min(grid::unit.c(width, height)) * 0.5,
        gp = grid::gpar(fill = nes_colors(nes_matrix[i, j]), col = "white", lwd = 0.5)
      )
    },
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 9),
    column_names_gp = grid::gpar(fontsize = 9),
    column_title = paste(title_prefix, "- Focused Pathway Analysis"),
    left_annotation = ComplexHeatmap::rowAnnotation(
      Category = plot_data$category[match(rownames(nes_matrix), plot_data$pathway_clean)] %>% factor(),
      col = list(Category = category_colors),
      annotation_name_gp = grid::gpar(fontsize = 10),
      show_annotation_name = TRUE,
      annotation_width = grid::unit(1, "cm")
    ),
    heatmap_legend_param = list(
      title = "NES",
      title_gp = grid::gpar(fontsize = 10),
      labels_gp = grid::gpar(fontsize = 9),
      legend_height = grid::unit(4, "cm")
    )
  )
  ComplexHeatmap::draw(ht, annotation_legend_list = list(size_legend))
  dev.off()
  png(paste0(output_path, "_focused_heatmap.png"), width = 14, height = 10, units = "in", res = 300)
  ComplexHeatmap::draw(ht, annotation_legend_list = list(size_legend))
  dev.off()
  return(ht)
}

create_focused_pathway_heatmap_nes_based <- function(focused_data, title_prefix, output_path) {
  if (nrow(focused_data) == 0) return(NULL)
  if ("p.adjust" %in% colnames(focused_data) && !"padj" %in% colnames(focused_data)) {
    focused_data <- focused_data %>% dplyr::rename(padj = p.adjust)
  }
  focused_pathways_order <- readr::read_csv(cfg$focused_pathways_path, col_names = c("Category", "Pathway"))
  focused_pathways_order <- focused_pathways_order %>% 
    dplyr::filter(!is.na(Category) & !is.na(Pathway) & Category != "" & Pathway != "") %>%
    dplyr::mutate(
      MSigDB_Name = paste0("HALLMARK_", sapply(Pathway, function(x) {
        x %>%
          stringr::str_to_upper() %>%
          stringr::str_replace_all("[^A-Z0-9]", "_") %>%
          stringr::str_replace_all("_+", "_") %>%
          stringr::str_remove("^_|_$")
      })),
      pathway_clean = stringr::str_to_title(stringr::str_replace_all(Pathway, "_", " "))
    )
  plot_data <- focused_data %>%
    dplyr::select(Description, NES, padj, cell_type, method, category) %>%
    dplyr::mutate(
      comparison = paste(cell_type, method, sep = "_"),
      abs_NES = abs(NES),
      color_category = dplyr::case_when(
        padj > 0.05 ~ "Non-significant",
        NES > 0 & padj <= 0.05 ~ "Significant Positive",
        NES < 0 & padj <= 0.05 ~ "Significant Negative",
        TRUE ~ "Non-significant"
      ),
      pathway_clean = stringr::str_remove(Description, "HALLMARK_"),
      pathway_clean = stringr::str_replace_all(pathway_clean, "_", " "),
      pathway_clean = stringr::str_to_title(pathway_clean)
    ) %>%
    dplyr::left_join(
      focused_pathways_order %>% dplyr::select(MSigDB_Name, Category) %>% dplyr::rename(Description = MSigDB_Name, category_ordered = Category),
      by = "Description"
    ) %>%
    dplyr::mutate(category = ifelse(!is.na(category_ordered), category_ordered, category)) %>%
    dplyr::select(-category_ordered)
  pathway_order <- focused_pathways_order %>% dplyr::filter(pathway_clean %in% plot_data$pathway_clean) %>% dplyr::pull(pathway_clean)
  additional_pathways <- setdiff(unique(plot_data$pathway_clean), pathway_order)
  pathway_order <- c(pathway_order, additional_pathways)
  nes_matrix <- plot_data %>% dplyr::select(pathway_clean, comparison, NES) %>% tidyr::pivot_wider(names_from = comparison, values_from = NES, values_fill = 0) %>% tibble::column_to_rownames("pathway_clean") %>% as.matrix()
  abs_nes_matrix <- plot_data %>% dplyr::select(pathway_clean, comparison, abs_NES) %>% tidyr::pivot_wider(names_from = comparison, values_from = abs_NES, values_fill = 0) %>% tibble::column_to_rownames("pathway_clean") %>% as.matrix()
  significance_matrix <- plot_data %>% dplyr::select(pathway_clean, comparison, color_category) %>% tidyr::pivot_wider(names_from = comparison, values_from = color_category, values_fill = "Non-significant") %>% tibble::column_to_rownames("pathway_clean") %>% as.matrix()
  common_pathways <- Reduce(intersect, list(rownames(nes_matrix), rownames(abs_nes_matrix), rownames(significance_matrix)))
  common_comparisons <- Reduce(intersect, list(colnames(nes_matrix), colnames(abs_nes_matrix), colnames(significance_matrix)))
  nes_matrix <- nes_matrix[common_pathways, common_comparisons, drop = FALSE]
  abs_nes_matrix <- abs_nes_matrix[common_pathways, common_comparisons, drop = FALSE]
  significance_matrix <- significance_matrix[common_pathways, common_comparisons, drop = FALSE]
  if (nrow(nes_matrix) == 0 || ncol(nes_matrix) == 0) return(NULL)
  significance_colors <- c("Significant Positive" = "red", "Significant Negative" = "blue", "Non-significant" = "gray")
  max_abs_nes <- max(abs_nes_matrix, na.rm = TRUE)
  size_fun <- function(x) { pmax(0.1, pmin(1, x / max_abs_nes)) }
  size_legend_values <- c(0.5, 1.0, 1.5, 2.0)
  size_legend_values <- size_legend_values[size_legend_values <= max_abs_nes]
  if (length(size_legend_values) == 0) size_legend_values <- c(max_abs_nes/2, max_abs_nes)
  size_legend <- ComplexHeatmap::Legend(title = "|NES|", title_gp = grid::gpar(fontsize = 10), labels = as.character(round(size_legend_values, 2)), labels_gp = grid::gpar(fontsize = 9), type = "points", pch = 16, size = grid::unit(size_fun(size_legend_values) * 8, "mm"), legend_height = grid::unit(3, "cm"))
  color_legend <- ComplexHeatmap::Legend(title = "Enrichment", title_gp = grid::gpar(fontsize = 10), labels = c("NES>0 & padj<0.05", "NES<0 & padj<0.05", "padj≥0.05"), labels_gp = grid::gpar(fontsize = 9), legend_gp = grid::gpar(fill = c("red", "blue", "gray")), type = "points", pch = 16, size = grid::unit(4, "mm"), legend_height = grid::unit(3, "cm"))
  pdf(paste0(output_path, "_nes_based_heatmap.pdf"), width = 14, height = 10)
  ht <- ComplexHeatmap::Heatmap(
    nes_matrix,
    name = "NES_dummy",
    show_heatmap_legend = FALSE,
    rect_gp = grid::gpar(type = "none"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      dot_color <- significance_colors[ significance_matrix[i, j] ]
      dot_size <- size_fun(abs_nes_matrix[i, j])
      grid::grid.circle(x = x, y = y, r = dot_size * min(grid::unit.c(width, height)) * 0.4, gp = grid::gpar(fill = dot_color, col = "white", lwd = 0.5))
    },
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 9),
    column_names_gp = grid::gpar(fontsize = 9),
    column_title = paste(title_prefix, "- NES-based Focused Pathway Analysis")
  )
  ComplexHeatmap::draw(ht, annotation_legend_list = list(size_legend, color_legend), annotation_legend_side = "right", merge_legend = FALSE)
  dev.off()
  png(paste0(output_path, "_nes_based_heatmap.png"), width = 14, height = 10, units = "in", res = 300)
  ComplexHeatmap::draw(ht, annotation_legend_list = list(size_legend, color_legend), annotation_legend_side = "right", merge_legend = FALSE)
  dev.off()
  return(ht)
}

# ================================================================================
# 7. Generate Visualizations ####

cat("\n=== Generating Visualizations ===\n")

# GO Visualizations - Separate by cell type
cat("Creating", GO_DISPLAY_NAME, "visualizations for each cell type...\n")
if (!is.null(combined_go_results) && nrow(combined_go_results) > 0) {
  go_cell_types <- unique(combined_go_results$cell_type)
  go_methods <- unique(combined_go_results$method)
  cat("Creating", GO_DISPLAY_NAME, "plots for", length(go_cell_types), "cell types and", length(go_methods), "methods\n")
  for (cell_type in go_cell_types) {
    for (method in go_methods) {
      cell_method_data <- combined_go_results %>% dplyr::filter(cell_type == !!cell_type & method == !!method)
      if (nrow(cell_method_data) > 0) {
        clean_cell_name <- gsub("[^A-Za-z0-9]", "_", cell_type)
        output_prefix <- paste0("step4_pathway_enrichment/figures/", GO_DIR_NAME, "/", clean_cell_name, "_", method, "_", GO_DIR_NAME)
        cat("  Creating", GO_DISPLAY_NAME, "plots for", cell_type, "-", method, "(", nrow(cell_method_data), "pathways)\n")
        create_lollipop_plot_by_nes(cell_method_data, paste(GO_DISPLAY_NAME, ":", cell_type, "-", method), output_prefix)
        create_lollipop_plot_by_significance(cell_method_data, paste(GO_DISPLAY_NAME, ":", cell_type, "-", method), output_prefix)
      }
    }
  }
} else {
  cat("No", GO_DISPLAY_NAME, "results to visualize\n")
}

# KEGG Visualizations - Separate by cell type
cat("Creating KEGG visualizations for each cell type...\n")
if (!is.null(combined_kegg_results) && nrow(combined_kegg_results) > 0) {
  kegg_cell_types <- unique(combined_kegg_results$cell_type)
  kegg_methods <- unique(combined_kegg_results$method)
  cat("Creating KEGG plots for", length(kegg_cell_types), "cell types and", length(kegg_methods), "methods\n")
  for (cell_type in kegg_cell_types) {
    for (method in kegg_methods) {
      cell_method_data <- combined_kegg_results %>% dplyr::filter(cell_type == !!cell_type & method == !!method)
      if (nrow(cell_method_data) > 0) {
        clean_cell_name <- gsub("[^A-Za-z0-9]", "_", cell_type)
        output_prefix <- paste0("step4_pathway_enrichment/figures/KEGG/", clean_cell_name, "_", method, "_KEGG")
        cat("  Creating KEGG plots for", cell_type, "-", method, "(", nrow(cell_method_data), "pathways)\n")
        create_lollipop_plot_by_nes(cell_method_data, paste("KEGG:", cell_type, "-", method), output_prefix)
        create_lollipop_plot_by_significance(cell_method_data, paste("KEGG:", cell_type, "-", method), output_prefix)
      }
    }
  }
} else {
  cat("No KEGG results to visualize\n")
}

# Hallmark Visualizations - Separate by cell type
cat("Creating Hallmark visualizations for each cell type...\n")
if (!is.null(combined_hallmark_results) && nrow(combined_hallmark_results) > 0) {
  hallmark_cell_types <- unique(combined_hallmark_results$cell_type)
  hallmark_methods <- unique(combined_hallmark_results$method)
  cat("Creating Hallmark plots for", length(hallmark_cell_types), "cell types and", length(hallmark_methods), "methods\n")
  for (cell_type in hallmark_cell_types) {
    for (method in hallmark_methods) {
      cell_method_data <- combined_hallmark_results %>% dplyr::filter(cell_type == !!cell_type & method == !!method)
      if (nrow(cell_method_data) > 0) {
        clean_cell_name <- gsub("[^A-Za-z0-9]", "_", cell_type)
        output_prefix <- paste0("step4_pathway_enrichment/figures/Hallmark/", clean_cell_name, "_", method, "_Hallmark")
        cat("  Creating Hallmark plots for", cell_type, "-", method, "(", nrow(cell_method_data), "pathways)\n")
        create_lollipop_plot_by_nes(cell_method_data, paste("Hallmark:", cell_type, "-", method), output_prefix)
        create_lollipop_plot_by_significance(cell_method_data, paste("Hallmark:", cell_type, "-", method), output_prefix)
      }
    }
  }
} else {
  cat("No Hallmark results to visualize\n")
}

# Focused Pathway Visualizations - Combined heatmap approach
cat("Creating focused pathway visualizations...\n")
if (!is.null(combined_focused_results) && nrow(combined_focused_results) > 0) {
  focused_cell_types <- unique(combined_focused_results$cell_type)
  focused_methods <- unique(combined_focused_results$method)
  cat("Creating focused pathway heatmaps for", length(focused_cell_types), "cell types and", length(focused_methods), "methods\n")
  create_focused_pathway_heatmap(combined_focused_results, "MSigDB Hallmark Pathways", "step4_pathway_enrichment/figures/Focused_Pathways/combined_focused")
  create_focused_pathway_heatmap_nes_based(combined_focused_results, "MSigDB Hallmark Pathways", "step4_pathway_enrichment/figures/Focused_Pathways/combined_focused")
  for (method in focused_methods) {
    method_data <- combined_focused_results %>% dplyr::filter(method == !!method)
    if (nrow(method_data) > 0) {
      create_focused_pathway_heatmap(method_data, paste("Focused Pathways", method), paste0("step4_pathway_enrichment/figures/Focused_Pathways/", method, "_focused"))
      create_focused_pathway_heatmap_nes_based(method_data, paste("Focused Pathways", method), paste0("step4_pathway_enrichment/figures/Focused_Pathways/", method, "_focused"))
    }
  }
} else {
  cat("No focused pathway results to visualize\n")
}

# ================================================================================
# 8. Comparative Analysis ####

cat("\n=== Performing Comparative Analysis ===\n")

compare_methods <- function(results_data, analysis_type) {
  if (nrow(results_data) == 0) return(NULL)
  if ("padj" %in% colnames(results_data) && !"p.adjust" %in% colnames(results_data)) {
    results_data <- results_data %>% dplyr::rename(p.adjust = padj)
    cat("Renamed 'padj' to 'p.adjust' for", analysis_type, "analysis\n")
  }
  if (!"p.adjust" %in% colnames(results_data)) {
    cat("Warning: No p.adjust column found in", analysis_type, "results\n")
    cat("Available columns:", paste(colnames(results_data), collapse = ", "), "\n")
    return(NULL)
  }
  comparison_summary <- results_data %>%
    dplyr::group_by(method, cell_type) %>%
    dplyr::summarise(
      total_pathways = dplyr::n(),
      significant_pathways = sum(p.adjust < 0.05, na.rm = TRUE),
      upregulated_pathways = sum(NES > 0 & p.adjust < 0.05, na.rm = TRUE),
      downregulated_pathways = sum(NES < 0 & p.adjust < 0.05, na.rm = TRUE),
      mean_nes = mean(abs(NES), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(analysis_type = analysis_type)
  return(comparison_summary)
}

go_comparison <- compare_methods(combined_go_results, GO_DIR_NAME)
kegg_comparison <- compare_methods(combined_kegg_results, "KEGG")
focused_comparison <- compare_methods(combined_focused_results, "Focused_MSigDB")
hallmark_comparison <- compare_methods(combined_hallmark_results, "Hallmark")

all_comparisons <- dplyr::bind_rows(go_comparison, kegg_comparison, focused_comparison, hallmark_comparison)
readr::write_csv(all_comparisons, "step4_pathway_enrichment/results/method_comparison_summary.csv")

if (nrow(all_comparisons) > 0) {
  comparison_plot <- ggplot2::ggplot(all_comparisons, ggplot2::aes(x = method, y = significant_pathways, fill = analysis_type)) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::facet_wrap(~analysis_type, scales = "free_y") +
    ggplot2::labs(title = "Comparison of Significant Pathways by Method and Analysis Type", x = "Method", y = "Number of Significant Pathways", fill = "Analysis Type") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"),
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9)
    )
  ggplot2::ggsave("step4_pathway_enrichment/figures/method_comparison.pdf", comparison_plot, width = 12, height = 8)
  ggplot2::ggsave("step4_pathway_enrichment/figures/method_comparison.png", comparison_plot, width = 12, height = 8, dpi = 300)
}

# ================================================================================
# 9. Session Info and Summary ####

session_info <- sessionInfo()
saveRDS(session_info, "step4_pathway_enrichment/results/session_info.rds")

analysis_summary <- data.frame(
  Analysis_Type = c(GO_DIR_NAME, "KEGG", "Focused_MSigDB", "Hallmark"),
  Total_Pathways = c(
    ifelse(is.null(combined_go_results), 0, nrow(combined_go_results)),
    ifelse(is.null(combined_kegg_results), 0, nrow(combined_kegg_results)),
    ifelse(is.null(combined_focused_results), 0, nrow(combined_focused_results)),
    ifelse(is.null(combined_hallmark_results), 0, nrow(combined_hallmark_results))
  ),
  Significant_Pathways = c(
    ifelse(is.null(combined_go_results), 0, sum(combined_go_results$p.adjust < 0.05, na.rm = TRUE)),
    ifelse(is.null(combined_kegg_results), 0, sum(combined_kegg_results$p.adjust < 0.05, na.rm = TRUE)),
    ifelse(is.null(combined_focused_results), 0, sum(ifelse("padj" %in% colnames(combined_focused_results), combined_focused_results$padj, combined_focused_results$p.adjust) < 0.05, na.rm = TRUE)),
    ifelse(is.null(combined_hallmark_results), 0, sum(ifelse("padj" %in% colnames(combined_hallmark_results), combined_hallmark_results$padj, combined_hallmark_results$p.adjust) < 0.05, na.rm = TRUE))
  ),
  Unique_Cell_Types = c(
    ifelse(is.null(combined_go_results), 0, length(unique(combined_go_results$cell_type))),
    ifelse(is.null(combined_kegg_results), 0, length(unique(combined_kegg_results$cell_type))),
    ifelse(is.null(combined_focused_results), 0, length(unique(combined_focused_results$cell_type))),
    ifelse(is.null(combined_hallmark_results), 0, length(unique(combined_hallmark_results$cell_type)))
  )
)

readr::write_csv(analysis_summary, "step4_pathway_enrichment/results/analysis_summary.csv")

cat("\n=== Analysis Complete ===\n")
cat("Results saved in: step4_pathway_enrichment/results/\n")
cat("Figures saved in: step4_pathway_enrichment/figures/\n")
cat("Summary:\n")
print(analysis_summary)

cat("\nScript completed successfully!\n")


