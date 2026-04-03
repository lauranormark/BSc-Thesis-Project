library(survival)
library(survminer)
library(dplyr)
library(ggplot2)


cell_types <- c("B_cells", "Cytotoxic_cells", "Dendritic_cells", "Mast_cells", "Neutrophils", "Eosinophils", "Macrophages", "NK_cells", "T_cells_CD4", "T_cells_CD8", "T_cells_gamma_delta", "T_regulatory_cells", "Macrophages_M1", "Macrophages_M2", "Endothelial", "Fibroblasts", "Monocytes", "Plasma_cells")  
# replace with all your actual cell type column names

# ── 1. EXTRACT COX RESULTS USING HIGH/LOW CUTPOINT ───────────────────────────
results <- list()

for (sex in c("Female", "Male")) {
  
  data_sex <- data_wide[data_wide$SEX == sex, ]
  
  for (cell in cell_types) {
    
    vars_needed <- c("OS_MONTHS", "OS_event", "AGE", cell)
    data_sub    <- data_sex[, vars_needed]
    data_sub    <- data_sub[complete.cases(data_sub), ]
    
    if (nrow(data_sub) < 10) next
    
    tryCatch({
      cut      <- surv_cutpoint(data_sub, time = "OS_MONTHS", event = "OS_event", variables = cell)
      data_cat <- surv_categorize(cut)
      data_cat$AGE     <- data_sub$AGE
      data_cat[[cell]] <- relevel(as.factor(data_cat[[cell]]), ref = "low")
      
      formula <- as.formula(paste("Surv(OS_MONTHS, OS_event) ~ AGE +", cell))
      cox_fit <- coxph(formula, data = data_cat)
      cox_sum <- summary(cox_fit)
      
      # ── USE paste0(cell, "high") TO MATCH THE ACTUAL ROW NAME ────────────
      row_name <- paste0(cell, "high")
      hr_row   <- cox_sum$conf.int[row_name, ]
      p_row    <- cox_sum$coefficients[row_name, ]
      
      results[[paste(cell, sex)]] <- data.frame(
        cell_type = cell,
        sex       = sex,
        HR        = as.numeric(hr_row["exp(coef)"]),
        lower     = as.numeric(hr_row["lower .95"]),
        upper     = as.numeric(hr_row["upper .95"]),
        p_value   = as.numeric(p_row["Pr(>|z|)"])
      )
    }, error = function(e) {
      message("Skipping ", cell, " (", sex, "): ", e$message)
    })
  }
}

# ── COMBINE AND ADD LABELS ─────────────────────────────────────────────────────
plot_df <- do.call(rbind, results)
rownames(plot_df) <- NULL

plot_df$HR      <- as.numeric(plot_df$HR)
plot_df$lower   <- as.numeric(plot_df$lower)
plot_df$upper   <- as.numeric(plot_df$upper)
plot_df$p_value <- as.numeric(plot_df$p_value)

plot_df$label <- paste0(
  "HR=",  round(plot_df$HR, 2),
  " (",   round(plot_df$lower, 2),
  "-",    round(plot_df$upper, 2), ")",
  "  p=", ifelse(plot_df$p_value < 0.001, "<0.001", round(plot_df$p_value, 3))
)

head(plot_df) # should now have values

# ── 3. PLOT FUNCTION ──────────────────────────────────────────────────────────
plot_forest <- function(df, sex_group, color) {
  ggplot(df, aes(x = HR, y = cell_type)) +
    geom_point(size = 3, color = color) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, color = color) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_text(aes(label = label), hjust = -0.05, size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.5))) +
    theme_classic() +
    labs(
      title = paste("Hazard Ratios by Cell Type in GBM —", sex_group, "(High vs Low Abundance)"),
      x     = "Hazard Ratio (95% CI)",
      y     = "Cell Type"
    ) +
    theme(axis.text.y = element_text(face = "bold", size = 10))
}


# ── 4. GENERATE AND SAVE PLOTS ────────────────────────────────────────────────

# First check the data looks correct
str(plot_df)
head(plot_df)

# Use filter() instead of bracket subsetting
p_female <- plot_forest(filter(plot_df, sex == "Female"), "Female", "#E64B35")
p_male   <- plot_forest(filter(plot_df, sex == "Male"),   "Male",   "#4DBBD5")

p_female
p_male

ggsave("Plots/forest_highlow_female.pdf", plot = p_female, width = 12, height = 7)
ggsave("Plots/forest_highlow_male.pdf",   plot = p_male,   width = 12, height = 7)