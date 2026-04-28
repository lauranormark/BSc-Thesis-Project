getwd()

library(survival)
library(survminer)
library(dplyr)
library(tidyr)
library(ggplot2)

# ── 1. LOAD & RESHAPE DATA ────────────────────────────────────────────────────

data_immunoabundance <- read.delim(
  "Data/merged_immunoabundance_clinical.txt",
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

cell_types_m <- c("Cytotoxic_cells",
                  "T_cells_CD4", "T_cells_CD8", "Fibroblasts", "NK_cells")

cell_order_m       <- rev(cell_types_m)
cell_types_clean_m <- gsub("_", ".", cell_types_m)
cell_order_clean_m <- gsub("_", ".", cell_order_m)

data_wide <- data_wide %>%
  rename_with(~ gsub("_", ".", .x), all_of(cell_types_m))

# ── 3. CATEGORIZE CELLS INTO HIGH/LOW ─────────────────────────────────────────

categorize_cells <- function(data_sex, cell_types_clean_m) {
  data_cat <- data_sex[, c("OS_MONTHS", "OS_event", "AGE", cell_types_clean_m)]
  data_cat <- data_cat[complete.cases(data_cat), ]
  
  for (cell in cell_types_clean_m) {
    cut <- surv_cutpoint(data_cat, time = "OS_MONTHS", event = "OS_event", variables = cell)
    categorized      <- surv_categorize(cut)
    data_cat[[cell]] <- relevel(as.factor(categorized[[cell]]), ref = "low")
  }
  return(data_cat)
}

data_male <- categorize_cells(
  data_wide[data_wide$SEX == "Male", ],
  cell_types_clean_m
)

# ── 4. CREATE AGE GROUP, DROP MISSING, REMOVE CONTINUOUS AGE ─────────────────

data_male <- data_male %>%
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

# ── 5. FIT MULTIVARIATE COX MODEL ─────────────────────────────────────────────

formula_mva <- as.formula(
  paste("Surv(OS_MONTHS, OS_event) ~ age_group +",
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
plot_df_m$variable <- gsub("high$", " (high)", plot_df_m$variable)

# Include age group levels in the display order
cell_order_display <- c(
  paste0(gsub("\\.", " ", cell_order_clean_m), " (high)"),
  "age_group45-55",
  "age_groupAbove 55"
)

# Clean up age group label names
plot_df_m$variable <- gsub("age_group", "Age: ", plot_df_m$variable)

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
  geom_text(aes(label = label), hjust = -0.05, size = 3, color = "black") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.55))) +
  theme_classic() +
  labs(
    title = "Multivariate Cox — Males (age-adjusted)",
    x     = "Hazard Ratio (95% CI)",
    y     = NULL
  ) +
  theme(axis.text.y = element_text(face = "bold", size = 10))

forest_male

# ─────────────── 9. SAVE ───────────────────────────────────────────────────────────────────

ggsave("Plots/forest_MVA_male.pdf", plot = forest_male, width = 10, height = 7)
