getwd()

library(survival)
library(survminer)
library(dplyr)
library(tidyr)
library(ggplot2)

# ── 1. LOAD & RESHAPE DATA ────────────────────────────────────────────────────

data_immunoabundance <- read.delim(
  "Datafiler/merged_immunoabundance_clinical.txt",
  header = TRUE,
  check.names = FALSE
)

data_wide <- data_immunoabundance %>%
  dplyr::select(PATIENT_ID, SEX, AGE, OS_STATUS, OS_MONTHS, CELL_TYPE, ABUNDANCE) %>%
  pivot_wider(
    names_from  = CELL_TYPE,
    values_from = ABUNDANCE
  )

data_wide$OS_event <- ifelse(
  grepl("1|DECEASED|Dead", data_wide$OS_STATUS, ignore.case = TRUE), 1, 0
)

# ── 2. DEFINE CELL TYPES ──────────────────────────────────────────────────────

cell_types_m <- c("B_cells", "Cytotoxic_cells", "Dendritic_cells", "Mast_cells",
                  "Neutrophils", "Eosinophils", "Macrophages", "NK_cells",
                  "T_cells_CD4", "T_cells_CD8", "T_cells_gamma_delta",
                  "T_regulatory_cells", "Macrophages_M1", "Macrophages_M2",
                  "Endothelial", "Fibroblasts", "Monocytes", "Plasma_cells")

cell_order_m       <- rev(cell_types_m)
cell_types_clean_m <- gsub("_", ".", cell_types_m)
cell_order_clean_m <- gsub("_", ".", cell_order_m)

data_wide <- data_wide %>%
  rename_with(~ gsub("_", ".", .x), all_of(cell_types_m))


data_male <- data_wide %>%
  filter(SEX == "Male") %>%
  dplyr::select(OS_MONTHS, OS_event, AGE, all_of(cell_types_clean_m)) %>%
  na.omit()

data_male[cell_types_clean_m] <- scale(data_male[cell_types_clean_m])
data_male$AGE <- scale(data_male$AGE)

# ── 5. FIT MULTIVARIATE COX MODEL ─────────────────────────────────────────────

formula_mva <- as.formula(
  paste("Surv(OS_MONTHS, OS_event) ~ AGE +",              
        paste(cell_types_clean_m, collapse = " + "))
)

cox_male <- coxph(formula_mva, data = data_male)

summary(cox_male)

# ── 6. EXTRACT RESULTS (cell types only, exclude age_group rows) ──────────────

extract_cox_results <- function(cox_fit, label) {
  s <- summary(cox_fit)
  df <- data.frame(
    variable  = rownames(s$conf.int),
    group     = label,
    HR        = s$conf.int[, "exp(coef)"],
    lower     = s$conf.int[, "lower .95"],
    upper     = s$conf.int[, "upper .95"],
    p_value   = s$coefficients[, "Pr(>|z|)"],
    row.names = NULL
  )
  # Remove only the NA row, keep the age_group rows for display
  df <- df[!grepl("NA", df$variable), ]
  return(df)
}

plot_df_m <- extract_cox_results(cox_male, "Male")

# ── 7. CLEAN LABELS ───────────────────────────────────────────────────────────

# ── Clean labels ──────────────────────────────────────────────────────────────
plot_df_m$variable <- gsub("\\.", " ", plot_df_m$variable)

plot_df_m$label <- paste0(
  "HR=", round(plot_df_m$HR, 2),
  " (", round(plot_df_m$lower, 2), "-", round(plot_df_m$upper, 2), ")",
  "  p=", ifelse(plot_df_m$p_value < 0.001, "<0.001", round(plot_df_m$p_value, 3))
)

# ── 8. PLOT ───────────────────────────────────────────────────────────────────

forest_male <- ggplot(plot_df_m, aes(x = HR, y = variable)) +
  geom_point(size = 3, color = "#4DBBD5") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_text(aes(label = label), x = 23, hjust = 0, size = 3, color = "black") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.55))) +
  theme_classic() +
  labs(
    title = "Multivariate Cox — Males (age adjusted)",
    x     = "Hazard Ratio (95% CI)",
    y     = NULL
  ) +
  theme(axis.text.y = element_text(face = "bold", size = 10))

forest_male

# ─────────────── 9. SAVE ───────────────────────────────────────────────────────────────────

ggsave("Plots/forest_MVA_male.pdf", plot = forest_male, width = 10, height = 7)
