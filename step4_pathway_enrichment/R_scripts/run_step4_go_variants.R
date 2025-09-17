# =============================================================================
# Runner: Execute step4_comprehensive_gsea_analysis.R for multiple GO ontologies
# =============================================================================
# Reads comma-separated tokens from config/options.txt key `go_ontologies`
# and runs the comprehensive script once per token (e.g., ALL, BP, MF, CC).
# =============================================================================

rm(list = ls())
if (getwd() != "/data/y1005/GBmap_Pipeline_v1.0") {
  setwd("/data/y1005/GBmap_Pipeline_v1.0")
}

source("config/common_utils.R")
cfg <- read_global_config()
set.seed(cfg$seed)

tokens <- cfg$go_ontologies
# Normalize to character vector split by comma if needed
if (length(tokens) == 1) {
  tokens <- unlist(strsplit(tokens, ",", fixed = TRUE))
}
tokens <- toupper(trimws(tokens))

if (length(tokens) == 0) tokens <- c("ALL")

valid <- c("ALL", "BP", "MF", "CC")
tokens <- tokens[tokens %in% valid]
if (length(tokens) == 0) tokens <- c("ALL")

cat("Running comprehensive GSEA for GO ontologies:", paste(tokens, collapse = ", "), "\n")

for (tok in tokens) {
  GO_ONTOLOGY <- tok
  cat("\n=== Executing step4_comprehensive_gsea_analysis.R with GO_ONTOLOGY =", GO_ONTOLOGY, "===\n")
  source("step4_pathway_enrichment/R_scripts/step4_comprehensive_gsea_analysis.R")
}

cat("\nAll requested GO variants completed.\n")


