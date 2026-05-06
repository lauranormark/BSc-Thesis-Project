#setwd("SoloDev")
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

cell_types_f <- c("B_cells", "Cytotoxic_cells", "Dendritic_cells", "Mast_cells",
                  "Neutrophils", "Eosinophils", "Macrophages", "NK_cells",
                  "T_cells_CD4", "T_cells_CD8", "T_cells_gamma_delta",
                  "T_regulatory_cells", "Macrophages_M1", "Macrophages_M2",
                  "Endothelial", "Fibroblasts", "Monocytes", "Plasma_cells")


cell_order_f       <- rev(cell_types_f)
cell_types_clean_f <- gsub("_", ".", cell_types_f)
cell_order_clean_f <- gsub("_", ".", cell_order_f)

data_wide <- data_wide %>%
  rename_with(~ gsub("_", ".", .x), all_of(cell_types_f))


data_female <- data_wide %>%
  filter(SEX == "Female") %>%
  dplyr::select(OS_MONTHS, OS_event, AGE, all_of(cell_types_clean_f)) %>%
  na.omit()


# ── 4. CREATE AGE GROUP, DROP MISSING, REMOVE CONTINUOUS AGE ─────────────────

data_female <- data_female %>%
  mutate(
    age_group = case_when(
      AGE < 45             ~ "Below 45",
      AGE >= 45 & AGE < 55 ~ "45-55",
      AGE >= 55            ~ "Above 55",
      TRUE                 ~ NA_character_
    ),
    age_group = factor(age_group, levels = c("Below 45", "45-55", "Above 55"))
  ) %>%
  filter(!is.na(age_group)) %>%
  dplyr::select(-AGE)

#variablerna kan ha olika skala, det här göt att HR blir effekten per SD
data_female[cell_types_clean_f] <- scale(data_female[cell_types_clean_f])

# ── 5. FIT MULTIVARIATE COX MODEL ─────────────────────────────────────────────

formula_mva <- as.formula(
  paste("Surv(OS_MONTHS, OS_event) ~ age_group +",
        paste(cell_types_clean_f, collapse = " + "))
)

cox_female <- coxph(formula_mva, data = data_female)

summary(cox_female)

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

plot_df_f <- extract_cox_results(cox_female, "Female")

# ── 7. CLEAN LABELS ───────────────────────────────────────────────────────────

# ── Clean labels ──────────────────────────────────────────────────────────────
plot_df_f$variable <- gsub("\\.", " ", plot_df_f$variable)


# Include age group levels in the display order
cell_order_display <- c(
  paste0(gsub("\\.", " ", cell_order_clean_f), " (high)"),
  "age_group45-55",
  "age_groupAbove 55"
)

# Clean up age group label names
plot_df_f$variable <- gsub("age_group", "Age: ", plot_df_f$variable)

plot_df_f$label <- paste0(
  "HR=", round(plot_df_f$HR, 2),
  " (", round(plot_df_f$lower, 2), "-", round(plot_df_f$upper, 2), ")",
  "  p=", ifelse(plot_df_f$p_value < 0.001, "<0.001", round(plot_df_f$p_value, 3))
)

# ── 8. PLOT ───────────────────────────────────────────────────────────────────

forest_female <- ggplot(plot_df_f, aes(x = HR, y = variable)) +
  geom_point(size = 3, color = "#E64B35") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_text(aes(label = label), x = 1850, hjust = 0, size = 3, color = "black") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.55))) +
  theme_classic() +
  labs(
    title = "Multivariate Cox — Females (age-adjusted)",
    x     = "Hazard Ratio (95% CI)",
    y     = NULL
  ) +
  theme(axis.text.y = element_text(face = "bold", size = 10))

forest_female

# ─────────────── 9. SAVE ───────────────────────────────────────────────────────────────────

ggsave("Plots/forest_MVA_female.pdf", plot = forest_female, width = 10, height = 7)
