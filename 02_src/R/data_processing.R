suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

obj <- readRDS("../../../Downloads/annotated_Sabah_data_21Oct2022.rds")

# check if normalized - yes
obj@assays$RNA@data
obj@assays$RNA@counts

#check if umap has been conducted - yes
obj@reductions

DimPlot(obj, reduction = "umap")

# ============================================================
# HostScope Immune Program Definitions
# ============================================================
# Each program represents a discrete immune behavior.
# Programs are non-overlapping by INTENT (not gene exclusivity).
# ============================================================

hostscope_programs <- list(
  
  IFN_I_Response = c(
    "ISG15", "IFI6", "IFIT1", "IFITM3", "MX1", "IRF7"
  ),
  
  Cytotoxic_Effector = c(
    "NKG7", "GNLY", "PRF1", "GZMB", "GZMA"
  ),
  
  Inflammatory_Response = c(
    "TNF", "IL1B", "CCL3", "CCL4"
  ),
  
  Immune_Checkpoint = c(
    "CTLA4", "PDCD1", "LAG3", "HAVCR2"
  ),
  
  CD4_Tr1_Regulatory = c(
    "IL10", "MAF", "ICOS", "CTLA4"
  ),
  
  Breg_Regulatory = c(
    "IL10", "CD24", "CD38", "MS4A1"
  ),
  
  Antigen_Presentation = c(
    "HLA-DRA", "HLA-DRB1", "HLA-DQB1", "CD74"
  ),
  
  Myeloid_Suppressive = c(
    "S100A8", "S100A9", "LILRB1", "CTSD"
  ),
  
  Fc_Opsonophagocytosis = c(
    "FCGR3A", "FCGR1A", "LST1", "TYROBP"
  ),
  
  Metabolic_Stress = c(
    "HIF1A", "LDHA", "BNIP3"
  ),
  
  Humoral_Bcell = c(
    "MS4A1", "CD79A", "IGHM", "IGKC"
  ),
  
  Platelet_QC = c(
    "PPBP", "PF4", "NRGN"
  )
)

# ============================================================
# Lineage-Restricted Program Map
# ============================================================
# Defines which immune programs are valid in each lineage.
# Prevents semantic leakage.
# ============================================================

lineage_program_map <- list(
  
  CD4_T = c(
    "IFN_I_Response",
    "CD4_Tr1_Regulatory",
    "Immune_Checkpoint"
  ),
  
  CD8_T = c(
    "IFN_I_Response",
    "Cytotoxic_Effector",
    "Immune_Checkpoint"
  ),
  
  NK = c(
    "IFN_I_Response",
    "Cytotoxic_Effector",
    "Immune_Checkpoint"
  ),
  
  GammaDelta_T = c(
    "Cytotoxic_Effector",
    "Antigen_Presentation"
  ),
  
  B = c(
    "Humoral_Bcell",
    "Breg_Regulatory",
    "Antigen_Presentation"
  ),
  
  Plasma = c(
    "Humoral_Bcell"
  ),
  
  Mono_CD14 = c(
    "Inflammatory_Response",
    "Myeloid_Suppressive",
    "Antigen_Presentation",
    "Fc_Opsonophagocytosis"
  ),
  
  Mono_CD16 = c(
    "Inflammatory_Response",
    "Myeloid_Suppressive",
    "Fc_Opsonophagocytosis"
  ),
  
  DC = c(
    "Antigen_Presentation",
    "IFN_I_Response"
  ),
  
  Platelet = c(
    "Platelet_QC"
  )
)

celltype_to_lineage <- c(
  "CD4 T cells"       = "CD4_T",
  "CD8 T cells"       = "CD8_T",
  "NK cells"          = "NK",
  "NKT cells"         = "NK",          # innate-like
  "γδ T cells"        = "GammaDelta_T",
  "B cells"           = "B",
  "Plasma cells"      = "Plasma",
  "CD14 Monocytes"    = "Mono_CD14",
  "CD16 Monocytes"    = "Mono_CD16",
  "cDCs"              = "DC",
  "pDCs"              = "DC",
  "Platelets"         = "Platelet",
  
  # Optional / QC-only
  "Proliferating cells" = "Proliferating",
  "HSPCs"               = "HSPC",
  "Unknown"             = "Unknown"
)

obj$lineage <- unname(celltype_to_lineage[obj@meta.data$celltype])

table(obj$lineage, useNA = "ifany")




# ------------------------------------------------------------
# Lineage-Restricted Program Scoring
# ------------------------------------------------------------


DefaultAssay(obj) <- "RNA"

for (lineage_name in names(lineage_program_map)) {
  
  message("Scoring lineage: ", lineage_name)
  
  lineage_cells <- rownames(obj@meta.data)[
    obj@meta.data$lineage == lineage_name
  ]
  
  
  if (length(lineage_cells) < 20) next
  
  valid_programs <- lineage_program_map[[lineage_name]]
  
  program_list <- lapply(
    hostscope_programs[valid_programs],
    function(g) intersect(g, rownames(obj))
  )
  
  keep <- lengths(program_list) >= 2
  program_list <- program_list[keep]
  program_names <- names(program_list)
  
  if (length(program_list) == 0) next
  
  obj <- AddModuleScore(
    object   = obj,
    features = program_list,
    name     = paste0("TMP_", lineage_name, "_"),
    assay    = "RNA",
    layer    = "data",      # 👈 THIS IS THE FIX
    cells    = lineage_cells
  )
  
  tmp_cols <- grep(
    paste0("^TMP_", lineage_name, "_"),
    colnames(obj@meta.data),
    value = TRUE
  )
  
  final_names <- paste0(lineage_name, "_", program_names)
  
  colnames(obj@meta.data)[match(tmp_cols, colnames(obj@meta.data))] <-
    final_names
}



# ------------------------------------------------------------
# Save scored object
# ------------------------------------------------------------
saveRDS(obj, "../../../Downloads/pbmc_program_scored.rds")
