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

# ── 2. DEFINE CELL TYPES & SHARED ORDER ───────────────────────────────────────

cell_types <- c("B_cells", "Cytotoxic_cells", "Dendritic_cells", "Mast_cells",
                "Neutrophils", "Eosinophils", "Macrophages", "NK_cells",
                "T_cells_CD4", "T_cells_CD8", "T_cells_gamma_delta",
                "T_regulatory_cells", "Macrophages_M1", "Macrophages_M2",
                "Endothelial", "Fibroblasts", "Monocytes", "Plasma_cells")

# Shared display order — must match script 07 exactly
cell_order <- rev(cell_types)

# Rename underscores to dots to avoid parsing issues
cell_types_clean <- gsub("_", ".", cell_types)
cell_order_clean <- gsub("_", ".", cell_order)

data_wide <- data_wide %>%
  rename_with(~ gsub("_", ".", .x), all_of(cell_types))

# ── 3. CATEGORIZE CELL TYPES INTO HIGH/LOW PER SEX ───────────────────────────

categorize_cells <- function(data_sex, cell_types_clean) {
  data_cat <- data_sex[, c("OS_MONTHS", "OS_event", "AGE", cell_types_clean)]
  data_cat  <- data_cat[complete.cases(data_cat), ]
  
  for (cell in cell_types_clean) {
    cut <- surv_cutpoint(data_cat, time = "OS_MONTHS", event = "OS_event", variables = cell)
    categorized      <- surv_categorize(cut)
    data_cat[[cell]] <- relevel(as.factor(categorized[[cell]]), ref = "low")
  }
  return(data_cat)
}

data_female <- categorize_cells(data_wide[data_wide$SEX == "Female", ], cell_types_clean)
data_male   <- categorize_cells(data_wide[data_wide$SEX == "Male", ],   cell_types_clean)

# ── 4. FIT COX MODELS ─────────────────────────────────────────────────────────

formula_all <- as.formula(
  paste("Surv(OS_MONTHS, OS_event) ~ AGE +",
        paste(cell_types_clean, collapse = " + "))
)

cox_female <- coxph(formula_all, data = data_female)
cox_male   <- coxph(formula_all, data = data_male)

summary(cox_female)
summary(cox_male)

# ── 5. EXTRACT RESULTS ────────────────────────────────────────────────────────

extract_cox_results <- function(cox_fit, sex_label) {
  s <- summary(cox_fit)
  data.frame(
    variable  = rownames(s$conf.int),
    sex       = sex_label,
    HR        = s$conf.int[, "exp(coef)"],
    lower     = s$conf.int[, "lower .95"],
    upper     = s$conf.int[, "upper .95"],
    p_value   = s$coefficients[, "Pr(>|z|)"],
    row.names = NULL
  )
}

results_female <- extract_cox_results(cox_female, "Female")
results_male   <- extract_cox_results(cox_male,   "Male")
plot_df        <- rbind(results_female, results_male)

# ── 6. CLEAN UP VARIABLE NAMES FOR DISPLAY ────────────────────────────────────

plot_df$variable <- gsub("\\.", " ", plot_df$variable)
plot_df$variable <- gsub("high$", " (high)", plot_df$variable)
plot_df$variable <- gsub("^AGE$", "Age (per 1-year increase)", plot_df$variable)

# Assign group for faceting
plot_df$group <- ifelse(
  grepl("Age", plot_df$variable), "Clinical", "Immune Cell"
)
plot_df$group <- factor(plot_df$group, levels = c("Clinical", "Immune Cell"))

# Apply shared cell order — must match cell_order from script 07
cell_order_display <- paste0(gsub("\\.", " ", cell_order_clean), " (high)")
plot_df$variable <- factor(
  plot_df$variable,
  levels = c("Age (per 1-year increase)", cell_order_display)
)

# Add HR labels
plot_df$label <- paste0(
  "HR=", round(plot_df$HR, 2),
  " (", round(plot_df$lower, 2), "-", round(plot_df$upper, 2), ")",
  "  p=", ifelse(plot_df$p_value < 0.001, "<0.001", round(plot_df$p_value, 3))
)

# ── 7. PLOT FUNCTION ──────────────────────────────────────────────────────────

plot_forest <- function(df, sex_group, color) {
  ggplot(df, aes(x = HR, y = variable, color = group)) +
    geom_point(size = 3) +
    geom_errorbar(aes(xmin = lower, xmax = upper),
                  height = 0.2, orientation = "y") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_text(aes(label = label), hjust = -0.05, size = 3, color = "black") +
    scale_color_manual(values = c("Clinical" = "#7E6148", "Immune Cell" = color)) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.55))) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    theme_classic() +
    labs(
      title = paste("Hazard Ratios in GBM —", sex_group, "(High vs Low Abundance)"),
      x     = "Hazard Ratio (95% CI)",
      y     = NULL,
      color = "Variable type"
    ) +
    theme(
      axis.text.y      = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text       = element_text(face = "bold")
    )
}

# ── 8. GENERATE AND SAVE ──────────────────────────────────────────────────────

p_female <- plot_forest(filter(plot_df, sex == "Female"), "Female", "#E64B35")
p_male   <- plot_forest(filter(plot_df, sex == "Male"),   "Male",   "#4DBBD5")

p_female
p_male

ggsave("Plots/forest_female.pdf", plot = p_female, width = 13, height = 8)
ggsave("Plots/forest_male.pdf",   plot = p_male,   width = 13, height = 8)
