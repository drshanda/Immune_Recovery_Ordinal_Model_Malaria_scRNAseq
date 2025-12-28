library(tibble)
library(pheatmap)
library(tidyverse)
library(broom)

obj <- readRDS("../../../Downloads/pbmc_program_scored.rds")

meta <- obj@meta.data %>%
  dplyr::select(ID, timepoint, celltype)

comp <- meta %>%
  count(ID, timepoint, celltype) %>%
  group_by(ID, timepoint) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

comp_patient <- comp %>%
  group_by(ID, celltype) %>%
  summarise(prop = mean(prop), .groups = "drop") %>%
  pivot_wider(names_from = celltype, values_from = prop, values_fill = 0)

mat <- as.matrix(comp_patient[, -1])
rownames(mat) <- comp_patient$ID


pheatmap(
  mat,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  color = colorRampPalette(c("navy", "white", "firebrick"))(100),
  main = "PBMC Composition per Patient",
  filename = "../../03_results/figures/composition/PBMC Composition per Patient.png"
)


df <- read.csv("../../03_results/tables/per_sample_shap_program_vs_composition.csv")

wilcox_test <- tidy(wilcox.test(
  df$mean_abs_shap_program,
  df$mean_abs_shap_composition,
  paired = TRUE,
  alternative = "greater"  # programs > composition
))

write_csv(wilcox_test, "../../03_results/tables/Wilcox_test_results.csv")

df_long <- df %>%
  pivot_longer(
    cols = c(mean_abs_shap_program, mean_abs_shap_composition),
    names_to = "feature_type",
    values_to = "mean_abs_shap"
  ) %>%
  mutate(
    feature_type = recode(
      feature_type,
      mean_abs_shap_program = "Immune Programs",
      mean_abs_shap_composition = "Cell Composition"
    )
  )


ggplot(df_long,
       aes(x = feature_type,
           y = mean_abs_shap,
           group = sample_id)) +
  geom_line(alpha = 0.3, color = "gray50") +
  geom_point(size = 2, alpha = 0.8) +
  geom_boxplot(
    aes(group = feature_type),
    outlier.shape = NA,
    alpha = 0.3,
    width = 0.5
  ) +
  labs(
    title = "Per-sample Model Reliance: Programs vs Composition",
    y = "Mean |SHAP|",
    x = NULL
  ) +
  theme_classic()
ggsave("../../03_results/figures/composition/Per-sample Model Reliance: Programs vs Composition.png")



