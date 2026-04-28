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

data_wide$AGE_GROUPS <- cut(data_wide$AGE,
                            breaks = c(0, 45, 56, Inf),
                            labels = c("<45", "45-55", "55<"),
                            right = FALSE)

data_wide$OS_event <- ifelse(
  grepl("1|DECEASED|Dead", data_wide$OS_STATUS, ignore.case = TRUE), 1, 0
)
# --------------FEMALES --------------------------------------------------------
# -------------SELEKTERA FEMALE CELLS ------------------------------------------
#Rensa denna baserat på relevanta celler
cell_types_f <- c("B_cells", "Cytotoxic_cells", "Dendritic_cells", "Mast_cells",
                  "Neutrophils", "Eosinophils", "Macrophages", "NK_cells",
                  "T_cells_CD4", "T_cells_CD8", "T_cells_gamma_delta",
                  "T_regulatory_cells", "Macrophages_M1", "Macrophages_M2",
                  "Endothelial", "Fibroblasts", "Monocytes")

# Shared display order — must match script 07 exactly
cell_order_f <- rev(cell_types_f)

# Rename underscores to dots to avoid parsing issues
cell_types_clean_f <- gsub("_", ".", cell_types_f)
cell_order_clean_f <- gsub("_", ".", cell_order_f)

data_wide <- data_wide %>%
  rename_with(~ gsub("_", ".", .x), all_of(cell_types_f))

# --------------KATEGORISERA CELLER I HIGH/LOW ---------------------------------
categorize_cells <- function(data_sex, cell_types_clean_f) {
  data_cat <- data_sex[, c("OS_MONTHS", "OS_event", "AGE_GROUPS", cell_types_clean_f)]
  data_cat  <- data_cat[complete.cases(data_cat), ]
  
  for (cell in cell_types_clean_f) {
    cut <- surv_cutpoint(data_cat, time = "OS_MONTHS", event = "OS_event", variables = cell)
    categorized      <- surv_categorize(cut)
    data_cat[[cell]] <- relevel(as.factor(categorized[[cell]]), ref = "low")
  }
  return(data_cat)
}

data_female_young  <- categorize_cells(
  data_wide[data_wide$SEX == "Female" & data_wide$AGE_GROUPS == "<45",  ],
  cell_types_clean_f
)

data_female_middle <- categorize_cells(
  data_wide[data_wide$SEX == "Female" & data_wide$AGE_GROUPS == "45-55", ],
  cell_types_clean_f
)

data_female_old    <- categorize_cells(
  data_wide[data_wide$SEX == "Female" & data_wide$AGE_GROUPS == "55<",  ],
  cell_types_clean_f
)

# --------------FITTING COX MODEL-----------------------------------------------
formula_all <- as.formula(
  paste("Surv(OS_MONTHS, OS_event) ~ +",
        paste(cell_types_clean_f, collapse = " + "))
)

# young group
cox_female_young <- coxph(formula_all, data = data_female_young)

# middle group
cox_female_middle <- coxph(formula_all, data = data_female_middle)

# old group
cox_female_old <- coxph(formula_all, data = data_female_old)


summary(cox_female_young)
summary(cox_female_middle)
summary(cox_female_old)


# --------------GETTING RESULTS FROM COX----------------------------------------

run_univariate_cox <- function(data_cat, cell_types, age_label) {
  results <- lapply(cell_types, function(cell) {
    formula <- as.formula(paste("Surv(OS_MONTHS, OS_event) ~", cell))
    fit <- tryCatch(coxph(formula, data = data_cat), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    extract_cox_results(fit, age_label)
  })
  do.call(rbind, Filter(Negate(is.null), results))
}

extract_cox_results <- function(cox_fit, age_label) {
  s <- summary(cox_fit)
  data.frame(
    variable  = rownames(s$conf.int),
    sex       = age_label,
    HR        = s$conf.int[, "exp(coef)"],
    lower     = s$conf.int[, "lower .95"],
    upper     = s$conf.int[, "upper .95"],
    p_value   = s$coefficients[, "Pr(>|z|)"],
    row.names = NULL
  )
}

plot_df_f_young  <- run_univariate_cox(data_female_young,  cell_types_clean_f, "Female <45")
plot_df_f_middle <- run_univariate_cox(data_female_middle, cell_types_clean_f, "Female 45-55")
plot_df_f_old    <- run_univariate_cox(data_female_old,    cell_types_clean_f, "Female 55<")

# --------------LABELS FOR EACH AGE GROUP --------------------------------------

apply_labels <- function(plot_df) {
  plot_df$variable <- gsub("\\.", " ", plot_df$variable)
  plot_df$variable <- gsub("high$", " (high)", plot_df$variable)
  
  cell_order_display_f <- paste0(gsub("\\.", " ", cell_order_clean_f), " (high)")
  plot_df$variable <- factor(plot_df$variable, levels = cell_order_display_f)
  
  plot_df$label <- paste0(
    "HR=", round(plot_df$HR, 2),
    " (", round(plot_df$lower, 2), "-", round(plot_df$upper, 2), ")",
    "  p=", ifelse(plot_df$p_value < 0.001, "<0.001", round(plot_df$p_value, 3))
  )
  return(plot_df)
}

plot_df_f_young  <- apply_labels(plot_df_f_young)
plot_df_f_middle <- apply_labels(plot_df_f_middle)
plot_df_f_old    <- apply_labels(plot_df_f_old)

# ------------PLOTTING----------------------------------------------------------

plot_forest <- function(df, age_group, color) {
  ggplot(df, aes(x = HR, y = variable)) +
    geom_point(size = 3, color = color) +
    geom_errorbarh(aes(xmin = lower, xmax = upper),   # horizontal error bars
                   height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_text(aes(label = label), hjust = -0.05, size = 3, color = "black") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.55))) +
    theme_classic() +
    labs(
      title = paste("Hazard Ratios in GBM —", age_group, "(High vs Low Abundance)"),
      x     = "Hazard Ratio (95% CI)",
      y     = NULL
    ) +
    theme(
      axis.text.y      = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text       = element_text(face = "bold")
    )
}

# Plot each age group separately
forest_young  <- plot_forest(plot_df_f_young,  "Females <45",    "#4DBBD5")
forest_middle <- plot_forest(plot_df_f_middle, "Females 45-55",  "#E64B35")
forest_old    <- plot_forest(plot_df_f_old,    "Females 55<",    "#3C5488")

forest_young
forest_middle
forest_old















# --------------MALES ----------------------------------------------------------

cell_types_m       <- c("B_cells", "Cytotoxic_cells", "Dendritic_cells", "Mast_cells",
                        "Neutrophils", "Eosinophils", "Macrophages", "NK_cells",
                        "T_cells_CD4", "T_cells_CD8", "T_cells_gamma_delta",
                        "T_regulatory_cells", "Macrophages_M1", "Macrophages_M2",
                        "Endothelial", "Fibroblasts", "Monocytes")

cell_order_m       <- rev(cell_types_m)
cell_types_clean_m <- gsub("_", ".", cell_types_m)
cell_order_clean_m <- gsub("_", ".", cell_order_m)

# --------------KATEGORISERA CELLER I HIGH/LOW ---------------------------------

data_male <- categorize_cells(
  data_wide[data_wide$SEX == "Male", ],
  cell_types_clean_m
)

# --------------UNIVARIATE COX -------------------------------------------------

plot_df_m <- run_univariate_cox(data_male, cell_types_clean_m, "Male")

# --------------LABELS ---------------------------------------------------------

apply_labels_m <- function(plot_df) {
  plot_df$variable <- gsub("\\.", " ", plot_df$variable)
  plot_df$variable <- gsub("high$", " (high)", plot_df$variable)
  
  cell_order_display_m <- paste0(gsub("\\.", " ", cell_order_clean_m), " (high)")
  plot_df$variable <- factor(plot_df$variable, levels = cell_order_display_m)
  
  plot_df$label <- paste0(
    "HR=", round(plot_df$HR, 2),
    " (", round(plot_df$lower, 2), "-", round(plot_df$upper, 2), ")",
    "  p=", ifelse(plot_df$p_value < 0.001, "<0.001", round(plot_df$p_value, 3))
  )
  return(plot_df)
}

plot_df_m <- apply_labels_m(plot_df_m)

# --------------PLOTTING -------------------------------------------------------

forest_male <- plot_forest(plot_df_m, "Males", "#4DBBD5")
forest_male

# --------------SAVE PLOTS -----------------------------------------------------

# Female age group forest plots
ggsave("Plots/forest_female_young.pdf",  plot = forest_young,  width = 10, height = 7)
ggsave("Plots/forest_female_middle.pdf", plot = forest_middle, width = 10, height = 7)
ggsave("Plots/forest_female_old.pdf",    plot = forest_old,    width = 10, height = 7)

# Male forest plot
ggsave("Plots/forest_male.pdf", plot = forest_male, width = 10, height = 7)