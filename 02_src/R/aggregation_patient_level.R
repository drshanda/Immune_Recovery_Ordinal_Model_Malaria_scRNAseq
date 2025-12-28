# ============================================================
# HostScope: Donor × Timepoint Aggregation
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)

# -----------------------------
# Inputs
# -----------------------------
# Load scored Seurat object
obj <- readRDS("../../../Downloads/pbmc_program_scored.rds")

# Required metadata columns
required_meta <- c("ID", "timepoint", "lineage")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing required metadata columns: ", paste(missing_meta, collapse = ", "))
}

md <- obj@meta.data %>%
  tibble::rownames_to_column("cell_barcode")

# -----------------------------
# Identify program score columns
# -----------------------------
# We consider a "program score column" any column that begins with one of your lineage names + "_"
# Example: "CD4_T_IFN_I_Response"
lineage_levels <- sort(unique(md$lineage))
lineage_prefix_pattern <- paste0("^(", paste(stringr::str_replace_all(lineage_levels, "([\\^\\$\\.|\\(\\)\\[\\]\\*\\+\\?\\\\])", "\\\\\\1"), collapse = "|"), ")_")

program_cols <- colnames(md)[str_detect(colnames(md), lineage_prefix_pattern)]

if (length(program_cols) == 0) {
  stop("No program score columns detected. Check that program scoring ran and columns are lineage-prefixed.")
}

# Optional: exclude QC-only platelet program from modeling features (keep composition separately)
# If you want to keep Platelet_QC as a feature, set this to FALSE
exclude_platelet_qc_program <- TRUE

if (exclude_platelet_qc_program) {
  program_cols <- program_cols[!str_detect(program_cols, "^Platelet_")]
  if (length(program_cols) == 0) stop("After excluding Platelet_*, no program columns remain.")
}

# -----------------------------
# Helper: aggregate numeric columns safely
# -----------------------------
safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

safe_q <- function(x, q = 0.9) {
  if (all(is.na(x))) return(NA_real_)
  as.numeric(quantile(x, probs = q, na.rm = TRUE, names = FALSE))
}

# -----------------------------
# A) Aggregate program scores to donor × timepoint (means by default)
# -----------------------------
# Notes:
# - Because scores are lineage-restricted, most cells will have NA for programs outside their lineage.
# - Aggregation should therefore be done on the cells where the feature is defined (na.rm=TRUE).

agg_programs_mean <- md %>%
  group_by(ID, timepoint) %>%
  summarise(across(all_of(program_cols), safe_mean), .groups = "drop")

# Optional: add additional summary stats (median, 90th percentile) if you want
compute_extra_summaries <- FALSE

if (compute_extra_summaries) {
  
  agg_programs_median <- md %>%
    group_by(ID, timepoint) %>%
    summarise(across(all_of(program_cols), safe_median), .groups = "drop") %>%
    rename_with(~ paste0(.x, "__median"), all_of(program_cols))
  
  agg_programs_q90 <- md %>%
    group_by(ID, timepoint) %>%
    summarise(across(all_of(program_cols), ~ safe_q(.x, 0.9)), .groups = "drop") %>%
    rename_with(~ paste0(.x, "__q90"), all_of(program_cols))
  
  agg_programs <- agg_programs_mean %>%
    left_join(agg_programs_median, by = c("ID", "timepoint")) %>%
    left_join(agg_programs_q90, by = c("ID", "timepoint"))
  
} else {
  agg_programs <- agg_programs_mean
}

# -----------------------------
# B) Add lineage composition features (recommended)
# -----------------------------
# Composition features help separate "state change" vs "cell mix change".
# These are fractions of each lineage within donor × timepoint.

composition <- md %>%
  count(ID, timepoint, lineage, name = "n_cells") %>%
  group_by(ID, timepoint) %>%
  mutate(total_cells = sum(n_cells),
         frac = n_cells / total_cells) %>%
  ungroup() %>%
  select(ID, timepoint, lineage, frac) %>%
  pivot_wider(
    names_from = lineage,
    values_from = frac,
    names_prefix = "comp_",
    values_fill = 0
  )

# -----------------------------
# C) Merge programs + composition
# -----------------------------
features <- agg_programs %>%
  left_join(composition, by = c("ID", "timepoint"))

# -----------------------------
# D) Optional: infected-only filters for your two locked models
# -----------------------------
# Pairwise D0 vs D28 and Ordinal D0/D7/D28 (infected only).
infected_timepoints <- c("Day0", "Day7", "Day28")

features_infected <- features %>%
  filter(timepoint %in% infected_timepoints)

features_pairwise <- features_infected %>%
  filter(timepoint %in% c("Day0", "Day28")) %>%
  mutate(y_binary = ifelse(timepoint == "Day0", 1L, 0L))  # 1=acute, 0=recovered

features_ordinal <- features_infected %>%
  mutate(y_stage = case_when(
    timepoint == "Day0" ~ 0L,
    timepoint == "Day7" ~ 1L,
    timepoint == "Day28" ~ 2L,
    TRUE ~ NA_integer_
  ))

# -----------------------------
# E) Save outputs
# -----------------------------

write.csv(features, "../../01_data/processed/hostscope_features_all.csv", row.names = FALSE)
write.csv(features_infected, "../../01_data/processed/hostscope_features_infected.csv", row.names = FALSE)
write.csv(features_pairwise, "../../01_data/processed/hostscope_features_pairwise_D0_vs_D28.csv", row.names = FALSE)
write.csv(features_ordinal, "../../01_data/processed/hostscope_features_ordinal_D0_D7_D28.csv", row.names = FALSE)

message("Saved feature matrices to ../../01_data/processed/")
