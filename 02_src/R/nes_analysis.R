library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(tibble)
library(tidyverse)
library(pheatmap)

obj <- readRDS("../../../Downloads/pbmc_program_scored.rds")


base_dir <- "../../03_results/figures/functional_analysis_plots"
dir.create(base_dir, showWarnings = FALSE)

safe_name <- function(x) {
  x |>
    gsub(" ", "_", x = _) |>
    gsub("[^A-Za-z0-9_]", "", x = _)
}


Idents(obj) <- "celltype"

pseudobulk_by_lineage <- function(obj, lineage) {
  obj_sub <- subset(obj, idents = lineage)
  
  AggregateExpression(
    obj_sub,
    group.by = c("ID", "timepoint"),
    assays = "RNA",
    slot = "counts",
    return.seurat = TRUE
  )
}

run_de_and_rank <- function(pb_obj,
                            ident_1 = "Day0",
                            ident_2 = "Day28",
                            logfc.threshold = 0.25,
                            min.pct = 0.2) {
  
  Idents(pb_obj) <- "timepoint"
  
  de <- FindMarkers(
    pb_obj,
    ident.1 = ident_1,
    ident.2 = ident_2,
    test.use = "wilcox",
    logfc.threshold = logfc.threshold,
    min.pct = min.pct
  )
  
  # Create ranked vector for GSEA
  gene_ranks <- de$avg_log2FC
  names(gene_ranks) <- rownames(de)
  gene_ranks <- sort(gene_ranks, decreasing = TRUE)
  
  list(
    de = de,
    gene_ranks = gene_ranks
  )
}

run_gsea_go <- function(gene_ranks) {
  
  gseGO(
    geneList = gene_ranks,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
  )
}

lineages <- c(
  "CD8 T cells",
  "CD4 T cells",
  "NK cells",
  "Plasma cells",
  "Platelets",
  "CD14 Monocytes"
)

results <- list()

for (lin in lineages) {
  
  message("Processing lineage: ", lin)
  
  lin_safe <- safe_name(lin)
  lin_dir <- file.path(base_dir, lin_safe)
  dir.create(lin_dir, showWarnings = FALSE)
  
  # ---- Pseudobulk ----
  pb <- pseudobulk_by_lineage(obj, lin)
  
  # ---- DE + ranking ----
  de_out <- run_de_and_rank(pb)
  gene_ranks <- de_out$gene_ranks
  
  # ---- GSEA ----
  ego <- run_gsea_go(gene_ranks)
  
  results[[lin]] <- list(
    pseudobulk = pb,
    de = de_out$de,
    gene_ranks = gene_ranks,
    gsea = ego
  )
  
  # ---- Plot only if enrichment exists ----
  if (!is.null(ego) && nrow(ego@result) > 0) {
    
    # ---- Dotplot ----
    p_dot <- dotplot(
      ego,
      showCategory = 12
    ) +
      ggtitle(paste(lin, ": GO Biological Processes (Recovery)"))
    
    ggsave(
      filename = file.path(lin_dir, "dotplot_GO_BP.png"),
      plot = p_dot,
      width = 7,
      height = 5,
      dpi = 300
    )
    
    # ---- Ridgeplot ----
    p_ridge <- ridgeplot(
      ego,
      showCategory = 12
    ) +
      ggtitle(paste(lin, ": GSEA Ridge Plot"))
    
    ggsave(
      filename = file.path(lin_dir, "ridgeplot_GO_BP.png"),
      plot = p_ridge,
      width = 7,
      height = 6,
      dpi = 300
    )
    
    # ---- Top pathway GSEA plot ----
    top_pathway <- ego@result$Description[1]
    
    p_gsea <- gseaplot2(
      ego,
      geneSetID = 1,
      title = paste(lin, "-", top_pathway)
    )
    
    ggsave(
      filename = file.path(lin_dir, "gseaplot_top_pathway.png"),
      plot = p_gsea,
      width = 7,
      height = 5,
      dpi = 300
    )
    
  } else {
    message("  No significant enrichment for ", lin)
  }
}


extract_nes_long <- function(results, top_n = 15) {
  
  nes_list <- lapply(names(results), function(lin) {
    
    ego <- results[[lin]]$gsea
    
    if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
    
    ego@result %>%
      arrange(desc(abs(NES))) %>%
      slice_head(n = top_n) %>%
      dplyr::select(
        pathway = Description,
        NES
      ) %>%
      mutate(lineage = lin)
  })
  
  bind_rows(nes_list)
}

nes_long <- extract_nes_long(results, top_n = 15)

nes_mat <- nes_long %>%
  pivot_wider(
    names_from = lineage,
    values_from = NES
  ) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

nes_mat_filt <- nes_mat[rowSums(!is.na(nes_mat)) >= 2, ]

row_vars <- apply(nes_mat_filt, 1, function(x) var(x, na.rm = TRUE))
nes_mat_filt <- nes_mat_filt[row_vars > 0, ]

scale_rows_bulletproof <- function(mat) {
  
  scaled <- t(apply(mat, 1, function(x) {
    
    # Keep only finite values
    x_finite <- x[is.finite(x)]
    
    # If fewer than 2 finite values → return zeros
    if (length(x_finite) < 2) {
      return(rep(0, length(x)))
    }
    
    sd_x <- sd(x_finite)
    
    # If zero or undefined variance → return zeros
    if (!is.finite(sd_x) || sd_x == 0) {
      return(rep(0, length(x)))
    }
    
    # Scale safely
    (x - mean(x_finite)) / sd_x
  }))
  
  rownames(scaled) <- rownames(mat)
  colnames(scaled) <- colnames(mat)
  
  scaled
}


nes_mat_scaled <- scale_rows_bulletproof(nes_mat_filt)

any(!is.finite(nes_mat_scaled))

pheatmap(
  nes_mat_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick"))(100),
  na_col = "grey90",
  border_color = NA,
  main = "Functional Program Enrichment Across Immune Lineages\n(Day28 vs Day0)",
  fontsize_row = 9,
  fontsize_col = 10,
  filename = file.path(base_dir, "NES_heatmap_Day28_vs_Day0.png"),
  width = 9,
  height = 10
)



write.csv(
  nes_mat_filt,
  file = file.path("../../03_results/tables/NES_matrix_filtered.csv")
)
