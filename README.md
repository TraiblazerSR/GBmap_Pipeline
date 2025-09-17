# GBmap Pipeline v1.0

### Overview
GBmap is a reproducible, stepwise single-cell RNA-seq analysis pipeline centered on a configurable target gene. Cells are split into high vs low expression groups using either Gaussian Mixture Modeling (GMM) or median thresholds and analyzed across visualization, differential expression, pathway enrichment, cell–cell communication (CellChat), and pseudotime (monocle3). Terminology throughout uses "high" and "low" by design.

The pipeline is designed to work with the GBmap dataset, a comprehensive single-cell and spatial atlas of IDH-wildtype glioblastoma (Ruiz-Moreno et al., 2025).

- Project root: `/GBmap_Pipeline_v1.0`
- Global config: `config/options.txt`
- Shared helpers: `config/common_utils.R`
- Steps: `step0` … `step6`, each with `R_scripts/`, `processed/`, `results/`, `figures/`, and `logs/` as applicable

### Key features
- Config-driven target gene and cell types across all steps
- High/low grouping by GMM or median (globally, and CellChat-specific override)
- Robust gene ID resolution (SYMBOL ↔ ENSEMBL; alias support)
- Consistent outputs and auto-created folders (no manual mkdir needed)
- Comprehensive pathway enrichment (GO, KEGG, MSigDB Hallmark, Focused)
- Rich CellChat analysis across unsplit and split annotation levels
- Pseudotime analyses including comparative high/low modeling

### Data
The GBmap dataset can be downloaded from:

**Primary source (CELLxGENE):**
`https://cellxgene.cziscience.com/collections/999f2a15-3d7e-440b-96ae-2c806799c08c`

**Alternative download (Baidu Netdisk):**
- **core.rds**
  - **Link:** https://pan.baidu.com/s/1DHWJ1QDanUl0L1G6qc_NxA?pwd=5ws4
  - **Passcode:** 5ws4
- **extended.h5ad**
  - **Link:** https://pan.baidu.com/s/1v58pmvCRCAWJjarAsWxIBA?pwd=apdq
  - **Passcode:** apdq

**Note:** Download the data in h5ad format and convert to RDS format before use. Place all data files in the `raw_data_to_read/data_GBmap/` directory.

### Configure
Edit `config/options.txt`:
- Data
  - `core_data_path`: path to Seurat object (RDS)
  - `annotation_cross_table_step1_path`: CSV mapping for Step 1
  - `annotation_cross_table_cellchat_path`: CSV mapping for Step 5
  - `annotation_cross_table_monocle_path`: CSV mapping for Step 6
  - `focused_pathways_path`: CSV file with focused pathways for Step 4
- Target gene and cells
  - `target_gene_symbol`, `target_gene_ensembl`, `target_gene_aliases`
  - `target_cell_types`: comma-separated list (e.g., Endothelial,TAM-BDM,TAM-MG)
  - `step6_target_cell_types`: cell types for Step 6 analysis (may differ from general target_cell_types)
- Grouping methods
  - `group_methods`: GMM | median | both
  - `both_primary_method`: GMM | median
  - `step5_group_methods`, `step5_both_primary_method`
  - `step5_levels_to_run`: tokens among `unsplit`, specific cell types, `allsplit`
- Enrichment and thresholds
  - `go_ontologies`: BP, MF, CC, or all
  - `gsea_celltypes_scope`: all | focused
  - `expression_threshold`: numeric cutoff for “expressed”
- Gene ID priority and reproducibility
  - `gene_id_priority`: e.g., symbol,ensembl
  - `seed`, `workers`

### How to run (per step)
Run each script with R (paths relative to project root):
- Step 0: Data exploration
  - `Rscript step0_data_exploration/R_scripts/step0_data_exploration.R`
- Step 1: Add annotations from cross-table
  - `Rscript step1_add_annotations_from_crosstable/R_scripts/step1_add_annotations_from_crosstable.R`
- Step 2: Visualization (UMAPs, counts, target gene)
  - `Rscript step2_visualization/R_scripts/visualize_annotations_and_target_gene.R`
- Step 3: Differential analysis (target high vs low)
  - FindMarkers: `Rscript step3_differential_analysis/R_scripts/findmarkers_target_gene_analysis.R`
  - Pseudobulk+DESeq2: `Rscript step3_differential_analysis/R_scripts/pseudobulk_target_gene_analysis.R`
- Step 4: Pathway enrichment (GO/KEGG/Hallmark/Focused)
  - `Rscript step4_pathway_enrichment/R_scripts/step4_comprehensive_gsea_analysis.R`
- Step 5: Cell–cell communication (CellChat)
  - Substep 1 (levels + dual UMAPs): `Rscript step5_cell_communication/R_scripts/substep1_add_annotations_and_umaps.R`
  - Substep 2 (CellChat per level): `Rscript step5_cell_communication/R_scripts/substep2_cellchat_analysis_all_levels.R`
  - Substep 3 (selectK): `Rscript step5_cell_communication/R_scripts/substep3_subpart1_selectK.R`
  - Substep 3 (patterns): edit k-values in `substep3_subpart2_patterns.R`, then run:
    - `Rscript step5_cell_communication/R_scripts/substep3_subpart2_patterns.R`
  - Substep 4 (bubble plots): `Rscript step5_cell_communication/R_scripts/substep4_bubble_plots.R`
- Step 6: Pseudotime (monocle3)
  - Pseudotime analysis: `Rscript step6_pseudotime_target_gene/R_scripts/substep1_pseudotime_analysis.R`

### What each step does

#### **Step 0: Data Exploration**
*Initial data inspection and quality assessment*

This foundational step provides comprehensive overview of your single-cell dataset using **Seurat** (Hao et al., 2021). It examines data structure, cell counts, feature distributions, and generates summary statistics to ensure data quality before downstream analysis.

**Key outputs:**
- Loads `core.rds`; writes structure summaries to `step0_data_exploration/results/`
- Dataset overview, cell/gene counts, metadata inspection
- Quality metrics and basic statistics

#### **Step 1: Annotation Enhancement** 
*Hierarchical cell type annotation mapping*

Enhances existing cell annotations by adding new annotation levels through cross-table mapping. Uses **dplyr** (Wickham et al., 2023) for efficient data manipulation and annotation joining based on existing annotation levels.

**Key outputs:**
- Reads cross-table; maps highest existing `annotation_level_{n}` to `annotation_level_{n+1}`
- Saves updated object: `step1_add_annotations_from_crosstable/processed/core_data.rds`
- Annotation mapping reports and validation summaries

#### **Step 2: Data Visualization**
*Comprehensive visualization of cell types and target gene expression*

Creates publication-ready visualizations using **ggplot2** (Wickham, 2016) and **Seurat** plotting functions. Generates UMAP embeddings, expression plots, and cell type distributions to provide intuitive data overview.

**Key outputs:**
- UMAPs for key `annotation_level_*` using **Seurat** UMAP implementation
- Cell-count barplots; target gene Feature/Violin/Dot plots (all cells and targeted subsets)
- Expression % by cell type and annotation summaries
- Multi-panel publication-ready figures

#### **Step 3: Differential Gene Expression Analysis**
*Statistical identification of target gene-associated expression changes*

Performs robust differential expression analysis using two complementary approaches: **Seurat**'s FindMarkers and pseudobulk analysis with **DESeq2** (Love et al., 2014). Cells are stratified by target gene expression (high vs low) using Gaussian Mixture Modeling (**mclust**; Scrucca et al., 2016) or median thresholds. For glioma subtype classification (AC-like, NPC-like, MES-like, OPC-like), refer to Neftel et al. (2019). For GMM clustering methodology, see Liu et al. (2022). For single-cell pseudobulk differential expression analysis, see Squair et al. (2021).

**Key outputs:**
- FindMarkers mode: per cell type, high vs low (GMM/median), filters and writes per-type and combined CSVs; jjVolcano-style multi-panel figures
- Pseudobulk mode: builds pseudo-samples (**DESeq2**); writes per-type and combined CSVs; jjVolcano figures
- Statistical summaries and effect size estimates

#### **Step 4: Pathway Enrichment Analysis**
*Functional interpretation through pathway analysis*

Conducts comprehensive pathway analysis using **clusterProfiler** (Wu et al., 2021) and **fgsea** (Korotkevich et al., 2021) for Gene Set Enrichment Analysis (GSEA). Integrates multiple databases including GO, KEGG, and MSigDB Hallmark pathways.

**Key outputs:**
- GSEA for GO (BP/MF/CC/ALL), KEGG, MSigDB Hallmark, and Focused (from config-specified path)
- Lollipop plots (by |NES| and significance), focused heatmaps (dot-size or NES-based) using **ComplexHeatmap** (Gu et al., 2016)
- Method comparison summary and enrichment statistics

#### **Step 5: Cell-Cell Communication Analysis**
*Intercellular signaling network inference*

Analyzes cell-cell communication patterns using **CellChat** (Jin et al., 2021), a comprehensive framework for inferring and analyzing intercellular communication networks. Examines ligand-receptor interactions across different cell groupings and expression states.

**Key outputs:**
- Substep 1: Adds `annotation_unsplit`; creates split levels per configured targets and methods:
  - Specific levels: `annotation_{CellType}_{Method}`
  - Combined: `annotation_allsplit_{Method}`
  - Records list: `results/new_levels_created.txt` and `results/levels_by_method.csv`
  - Dual UMAPs: `figures/dual_umaps/{unsplit|GMM|median}/`
- Substep 2: For each level:
  - Creates **CellChat** object, runs inference, exports communications CSV
  - Saves network circle plots, pathway hierarchy/circle/heatmap, LR contributions, signaling role scatter/heatmap
- Substep 3: selectK and pattern analyses (outgoing/incoming); edit k in `substep3_subpart2_patterns.R`
- Substep 4: Bubble plots for interactions, pathways, and selected LR pairs

#### **Step 6: Pseudotime Trajectory Analysis**
*Temporal dynamics and developmental trajectory inference*

Reconstructs cellular trajectories and pseudotemporal ordering using **monocle3** (Cao et al., 2019). Identifies genes that vary along pseudotime and performs trajectory-based differential expression analysis to understand dynamic processes.

**Key outputs:**
- Builds CDS from Seurat; roots graph at leaf with lowest target expression
- Per configured cell types:
  - Pseudotime UMAPs, target gene expression plots, trajectory visualization
  - Target co-varying genes, correlation heatmaps using **pheatmap** (Kolde, 2019)
  - Graph tests and gene module discovery with **monocle3**

### Outputs (examples)
- `step{n}/*/results/`: CSVs, combined summaries, logs, session info
- `step{n}/*/figures/`: PDFs/PNGs (UMAPs, barplots, jjVolcano, lollipops, heatmaps, CellChat plots, pseudotime)
- `step5_cell_communication/results/{method}/{level}/`: CellChat objects/plots per level
- `step6_pseudotime_target_gene/results/{cell_type}/`: correlations, graph tests, modules, CDS objects

### Conventions and notes
- Annotation columns: `annotation_level_{n}`, `annotation_unsplit`, `annotation_{CellType}_{Method}`, `annotation_allsplit_{Method}`
- Group labels: “high” and “low” only (not positive/negative)
- Scripts auto-create folders; avoid manual folder creation
- Figures are saved with `ggsave` where applicable

### Dependencies (R)
- Core: Seurat, SeuratObject, dplyr, ggplot2, ggrepel, patchwork, Matrix, readr, tidyr, stringr, purrr, randomcoloR, ggsci
- Differential: DESeq2, edgeR, limma
- Enrichment/plots: clusterProfiler, org.Hs.eg.db, DOSE, enrichplot, msigdbr, fgsea, ComplexHeatmap, circlize, UpSetR, viridis, RColorBrewer
- CellChat: CellChat, future, ggalluvial (optional: reticulate for UMAP similarity; currently not required)
- Pseudotime: monocle3, pheatmap
- Grouping: mclust

### Data expectations
- `core_data_path` points to a Seurat object; Step 1 requires a valid cross-table keyed to the highest existing `annotation_level_*`
- Focused pathway list: specified in config `focused_pathways_path`
- All data files should be placed in `raw_data_to_read/data_GBmap/` directory

### Reproducibility
- Global `seed` and `workers` from `config/options.txt` are respected across steps

### References and Citations

#### **GBmap Dataset**
Ruiz-Moreno, C., Salas, S. M., Samuelsson, E., Minaeva, M., Ibarra, I., Grillo, M., Brandner, S., Roy, A., Forsberg-Nilsson, K., Kranendonk, M. E. G., Theis, F. J., Nilsson, M., & Stunnenberg, H. G. (2025). Charting the Single-Cell and Spatial Landscape of IDH-Wildtype Glioblastoma with GBmap. *Neuro-oncology*, noaf113. Advance online publication. https://doi.org/10.1093/neuonc/noaf113

#### **Core Analysis Packages**

**Single-cell analysis framework:**
- Hao, Y., Hao, S., Andersen-Nissen, E., Mauck, W. M., 3rd, Zheng, S., Butler, A., Lee, M. J., Wilk, A. J., Darby, C., Zager, M., Hoffman, P., Stoeckius, M., Papalexi, E., Mimitou, E. P., Jain, J., Srivastava, A., Stuart, T., Fleming, L. M., Yeung, B., Rogers, A. J., ... Satija, R. (2021). Integrated analysis of multimodal single-cell data. *Cell*, 184(13), 3573-3587.e29. https://doi.org/10.1016/j.cell.2021.04.048

**Data manipulation and visualization:**
- Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag New York. https://ggplot2.tidyverse.org
- Wickham, H., François, R., Henry, L., Müller, K., & Vaughan, D. (2023). dplyr: A Grammar of Data Manipulation. R package version 1.1.4. https://CRAN.R-project.org/package=dplyr

**Differential expression analysis:**
- Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8
- Squair, J. W., Gautier, M., Kathe, C., Anderson, M. A., James, N. D., Hutson, T. H., Hudelle, R., Qaiser, T., Matson, K. J. E., Barraud, Q., Levine, A. J., La Manno, G., Skinnider, M. A., & Courtine, G. (2021). Confronting false discoveries in single-cell differential expression. *Nature communications*, 12(1), 5692. https://doi.org/10.1038/s41467-021-25960-2

**Biological classifications and methods:**
- Neftel, C., Laffy, J., Filbin, M. G., Hara, T., Shore, M. E., Rahme, G. J., Richman, A. R., Silverbush, D., Shaw, M. L., Hebert, C. M., Dewitt, J., Gritsch, S., Perez, E. M., Gonzalez Castro, L. N., Lan, X., Druck, N., Rodman, C., Dionne, D., Kaplan, A., Bertalan, M. S., … Suvà, M. L. (2019). An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma. *Cell*, 178(4), 835–849.e21. https://doi.org/10.1016/j.cell.2019.06.024

**Clustering and mixture modeling:**
- Scrucca, L., Fop, M., Murphy, T. B., & Raftery, A. E. (2016). mclust 5: Clustering, Classification and Density Estimation Using Gaussian Finite Mixture Models. *The R Journal*, 8(1), 289-317. https://doi.org/10.32614/RJ-2016-021
- Liu, T. C., Kalugin, P. N., Wilding, J. L., & Bodmer, W. F. (2022). GMMchi: gene expression clustering using Gaussian mixture modeling. *BMC bioinformatics*, 23(1), 457. https://doi.org/10.1186/s12859-022-05006-0

**Pathway and functional analysis:**
- Wu, T., Hu, E., Xu, S., Chen, M., Guo, P., Dai, Z., Feng, T., Zhou, L., Tang, W., Zhan, L., Fu, X., Liu, S., Bo, X., & Yu, G. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *Innovation*, 2(3), 100141. https://doi.org/10.1016/j.xinn.2021.100141
- Korotkevich, G., Sukhov, V., Budin, N., Shpak, B., Artyomov, M. N., & Sergushichev, A. (2021). Fast gene set enrichment analysis. *bioRxiv*. https://doi.org/10.1101/060012

**Visualization and heatmaps:**
- Gu, Z., Eils, R., & Schlesner, M. (2016). Complex heatmaps reveal patterns and correlations in multidimensional genomic data. *Bioinformatics*, 32(18), 2847-2849. https://doi.org/10.1093/bioinformatics/btw313
- Kolde, R. (2019). pheatmap: Pretty Heatmaps. R package version 1.0.12. https://CRAN.R-project.org/package=pheatmap

**Cell-cell communication:**
- Jin, S., Guerrero-Juarez, C. F., Zhang, L., Chang, I., Ramos, R., Kuan, C. H., Myung, P., Plikus, M. V., & Nie, Q. (2021). Inference and analysis of cell-cell communication using CellChat. *Nature Communications*, 12(1), 1088. https://doi.org/10.1038/s41467-021-21246-9

**Pseudotime and trajectory analysis:**
- Cao, J., Spielmann, M., Qiu, X., Huang, X., Ibrahim, D. M., Hill, A. J., Zhang, F., Mundlos, S., Christiansen, L., Steemers, F. J., Trapnell, C., & Shendure, J. (2019). The single-cell transcriptional landscape of mammalian organogenesis. *Nature*, 566(7745), 496-502. https://doi.org/10.1038/s41586-019-0969-x

### How to Cite This Pipeline
If you use this pipeline in your research, please cite:

1. The GBmap dataset: Ruiz-Moreno et al. (2025)
2. The core analysis packages listed above, particularly Seurat for single-cell analysis
3. Any specific methods/packages used in your analysis (DESeq2, CellChat, monocle3, etc.)

### Special Thanks
I would like to express my sincere gratitude to:

**Professor Zhan Ren-Ya**, Honorary Director of the Department of Neurosurgery, The First Affiliated Hospital, Zhejiang University School of Medicine.

Senior colleagues from the department: **Jin Lin-Chun/Kim**, **Yan Min**, **Shen Yuan-Yuan**, **Zhang Chao**, **Zhu Yu**, and other mentors.

Fellow researcher **He Hai-Feng**, **MUNGUR Rajneesh** and all other members of **Professor Zhan's research team** and **Neurosurgery Basic Lab of Prof.Zhan/Dr.Jin**.

**All medical staff** in the **Neurosurgery wards at both Qingchun Campus and Zhijiang Campus** of The First Affiliated Hospital of Zhejiang University.
