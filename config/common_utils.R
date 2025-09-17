suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(Seurat)
  library(mclust)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# =============================================================================
# SECTION 1: CONFIG LOADING ####
# =============================================================================

read_global_config <- function(config_path = file.path(getwd(), "config", "options.txt")) {
  if (!file.exists(config_path)) {
    stop("Config file not found: ", config_path)
  }
  lines <- readr::read_lines(config_path)
  lines <- lines[!grepl("^\\s*#", lines)]
  lines <- lines[nzchar(trimws(lines))]

  parse_line <- function(x) {
    kv <- strsplit(x, "=", fixed = TRUE)[[1]]
    key <- trimws(kv[1])
    value <- trimws(paste(kv[-1], collapse = "="))
    return(c(key, value))
  }

  kv_pairs <- lapply(lines, parse_line)
  cfg <- as.list(setNames(lapply(kv_pairs, `[[`, 2), vapply(kv_pairs, `[[`, "", 1)))

  # Vectorize selected keys
  vectorize <- function(x) {
    if (is.null(x) || nchar(x) == 0) return(character(0))
    out <- unlist(strsplit(x, ","))
    trimws(out)
  }

  cfg$target_cell_types <- vectorize(cfg$target_cell_types)
  cfg$go_ontologies <- vectorize(cfg$go_ontologies)
  cfg$gene_id_priority <- vectorize(cfg$gene_id_priority)
  cfg$target_gene_aliases <- vectorize(cfg$target_gene_aliases)

  # Step 6 specific target cell types (may differ from general target_cell_types)
  if (!is.null(cfg$step6_target_cell_types)) {
    cfg$step6_target_cell_types <- vectorize(cfg$step6_target_cell_types)
  } else {
    cfg$step6_target_cell_types <- cfg$target_cell_types
  }

  # Coerce some numerics
  cfg$expression_threshold <- suppressWarnings(as.numeric(cfg$expression_threshold))
  if (is.na(cfg$expression_threshold)) cfg$expression_threshold <- 0
  cfg$workers <- suppressWarnings(as.integer(cfg$workers))
  if (is.na(cfg$workers)) cfg$workers <- 1
  cfg$seed <- suppressWarnings(as.integer(cfg$seed))
  if (is.na(cfg$seed)) cfg$seed <- 1234

  # Defaults for optional keys
  if (is.null(cfg$gsea_celltypes_scope) || !tolower(cfg$gsea_celltypes_scope) %in% c("all", "focused")) {
    cfg$gsea_celltypes_scope <- "all"
  } else {
    cfg$gsea_celltypes_scope <- tolower(cfg$gsea_celltypes_scope)
  }

  # Step 5 specific options
  if (is.null(cfg$step5_group_methods)) cfg$step5_group_methods <- cfg$group_methods
  if (is.null(cfg$step5_both_primary_method)) cfg$step5_both_primary_method <- cfg$both_primary_method
  cfg$step5_levels_to_run <- vectorize(cfg$step5_levels_to_run)
  if (length(cfg$step5_levels_to_run) == 0) cfg$step5_levels_to_run <- c("unsplit", "allsplit")

  return(cfg)
}

get_project_root <- function(cfg) {
  file.path("/data/y1005", cfg$project_name)
}

# =============================================================================
# SECTION 2: DIRECTORY HELPERS ####
# =============================================================================

ensure_directories <- function(paths) {
  for (p in paths) {
    if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  }
}

# =============================================================================
# SECTION 3: ANNOTATION LEVEL HELPERS ####
# =============================================================================

get_highest_annotation_level <- function(meta_df) {
  level_cols <- grep("^annotation_level_\\d+$", colnames(meta_df), value = TRUE)
  if (length(level_cols) == 0) return(list(col = NULL, num = 0))
  nums <- as.integer(sub("annotation_level_", "", level_cols))
  max_idx <- which.max(nums)
  list(col = level_cols[max_idx], num = nums[max_idx])
}

next_annotation_col_name <- function(meta_df) {
  info <- get_highest_annotation_level(meta_df)
  paste0("annotation_level_", info$num + 1)
}

# Pick source and target columns from a cross-table based on what exists in meta.data
pick_mapping_columns <- function(cross_df, meta_df) {
  # Prefer mapping from the highest annotation level present in meta
  meta_info <- get_highest_annotation_level(meta_df)
  candidate_from <- c(meta_info$col, paste0("annotation_level_", rev(0:10)))
  from_col <- candidate_from[candidate_from %in% colnames(cross_df)][1]
  if (is.na(from_col) || is.null(from_col)) {
    # Fallback to any annotation_level_* present in cross_df
    cross_annos <- grep("^annotation_level_\\d+$", colnames(cross_df), value = TRUE)
    from_col <- cross_annos[1]
  }
  if (is.na(from_col) || is.null(from_col)) stop("No suitable 'from' column found in cross table")

  # Prefer known target columns
  preferred_to <- c("annotation_custom", "annotation_cellchat", "annotation_TAMs", "annotation_target", "annotation_new")
  to_col <- preferred_to[preferred_to %in% colnames(cross_df)][1]
  if (is.na(to_col) || is.null(to_col)) {
    # Fallback to the first non-from column
    to_col <- setdiff(colnames(cross_df), from_col)[1]
  }
  if (is.na(to_col) || is.null(to_col)) stop("No suitable 'to' column found in cross table")

  list(from = from_col, to = to_col)
}

add_mapped_annotation <- function(seurat_obj, from_col, to_col, mapping_df, new_col_name) {
  mapping_df <- mapping_df %>% dplyr::select(dplyr::all_of(from_col), dplyr::all_of(to_col)) %>% dplyr::distinct()
  map_vec <- setNames(mapping_df[[to_col]], mapping_df[[from_col]])
  source_vals <- as.character(seurat_obj@meta.data[[from_col]])
  seurat_obj@meta.data[[new_col_name]] <- unname(map_vec[source_vals])
  seurat_obj@meta.data[[new_col_name]] <- factor(seurat_obj@meta.data[[new_col_name]])
  seurat_obj
}

# =============================================================================
# SECTION 4: GENE ID AND EXPRESSION HELPERS ####
# =============================================================================

# Returns a feature name present in rownames(core_data) based on config priorities
identify_target_gene_feature <- function(seurat_obj, cfg) {
  features <- rownames(seurat_obj)
  # Try exact IDs from config
  if (!is.null(cfg$target_gene_ensembl) && cfg$target_gene_ensembl %in% features) return(cfg$target_gene_ensembl)
  if (!is.null(cfg$target_gene_symbol) && cfg$target_gene_symbol %in% features) return(cfg$target_gene_symbol)

  # Try aliases directly
  if (length(cfg$target_gene_aliases) > 0) {
    for (alias in cfg$target_gene_aliases) {
      if (alias %in% features) return(alias)
    }
  }

  # Try mapping symbols to Ensembl or vice versa using clusterProfiler::bitr
  try({
    # Heuristic: if features look like ENSG, map symbol->ENSEMBL; else map ENSEMBL->SYMBOL
    looks_ensembl <- mean(grepl("^ENSG", head(features, 100))) > 0.5
    if (!is.null(cfg$target_gene_symbol)) {
      if (looks_ensembl) {
        map_df <- clusterProfiler::bitr(cfg$target_gene_symbol, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Hs.eg.db)
        if (nrow(map_df) > 0) {
          ens <- unique(map_df$ENSEMBL)
          hit <- ens[ens %in% features]
          if (length(hit) > 0) return(hit[1])
        }
      }
    }
    if (!is.null(cfg$target_gene_ensembl)) {
      if (!looks_ensembl) {
        map_df <- clusterProfiler::bitr(cfg$target_gene_ensembl, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
        if (nrow(map_df) > 0) {
          sym <- unique(map_df$SYMBOL)
          hit <- sym[sym %in% features]
          if (length(hit) > 0) return(hit[1])
        }
      }
    }
  }, silent = TRUE)

  stop("Target gene not found in feature names using provided identifiers and aliases.")
}

get_expression_vector <- function(seurat_obj, feature) {
  expr <- tryCatch({
    Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "data")[feature, ]
  }, error = function(e) {
    Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts")[feature, ]
  })
  as.numeric(expr) %>% stats::setNames(colnames(seurat_obj))
}

# =============================================================================
# SECTION 5: GROUPING METHODS (GMM / MEDIAN) ####
# =============================================================================

split_by_gmm <- function(expr_vec) {
  # expr_vec: named numeric vector
  fit <- NULL
  model_types <- c("V", "E", "VVV", "EEE")
  for (m in model_types) {
    fit <- tryCatch({ mclust::Mclust(expr_vec, G = 2, modelNames = m) }, error = function(e) NULL)
    if (!is.null(fit) && !is.null(fit$G) && fit$G == 2) break
  }
  if (is.null(fit) || is.null(fit$G) || fit$G != 2) {
    # Fallback to median
    med <- stats::median(expr_vec)
    return(list(high = names(expr_vec)[expr_vec > med], low = names(expr_vec)[expr_vec <= med], threshold = med, method = "median_fallback"))
  }
  means <- fit$parameters$mean
  high_comp <- which.max(means)
  post <- fit$z
  high <- names(expr_vec)[post[, high_comp] > 0.5]
  low  <- names(expr_vec)[post[, high_comp] <= 0.5]
  list(high = high, low = low, threshold = NA_real_, method = "gmm", model = fit$modelName, means = means)
}

split_by_median <- function(expr_vec) {
  med <- stats::median(expr_vec)
  list(high = names(expr_vec)[expr_vec > med], low = names(expr_vec)[expr_vec <= med], threshold = med, method = "median")
}

# Create labels like "Endothelial_H2AJ_high" or "Endothelial_H2AJ_low"
label_high_low <- function(cell_barcodes, high_set, low_set, cell_type_name, gene_symbol) {
  labels <- rep(NA_character_, length(cell_barcodes))
  names(labels) <- cell_barcodes
  labels[names(labels) %in% high_set] <- paste0(cell_type_name, "_", gene_symbol, "_high")
  labels[names(labels) %in% low_set]  <- paste0(cell_type_name, "_", gene_symbol, "_low")
  labels
}

# Compute percent expressed (> threshold) by a grouping column
compute_expression_percentage <- function(seurat_obj, feature, group_col, threshold = 0) {
  expr_vec <- get_expression_vector(seurat_obj, feature)
  df <- seurat_obj@meta.data %>% dplyr::mutate(.group = .data[[group_col]], .expr = expr_vec[colnames(seurat_obj)])
  df %>% dplyr::filter(!is.na(.group)) %>% dplyr::group_by(.group) %>% dplyr::summarise(
    cells = dplyr::n(),
    positive = sum(.expr > threshold),
    percentage = round(positive / cells * 100, 2),
    .groups = 'drop'
  ) %>% dplyr::arrange(dplyr::desc(percentage))
}
