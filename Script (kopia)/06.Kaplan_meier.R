getwd()
#setwd("AgnesVärld")

#Install packages if needed:
install.packages("survival")
install.packages("survminer")

library(survival)
library(survminer)
library(dplyr)

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────

data_immunoabundance <- read.delim(
  "Filis/Datafiler/merged_immunoabundance_clinical.txt",
  header = TRUE,
  check.names = FALSE
)

# Kaplan-Meier behöver en rad per patient (inte long-format).
# Vi tar ut unika patienter med klinisk data.
clinical_unique <- data_immunoabundance %>%
  select(PATIENT_ID, AGE, SEX, OS_STATUS, OS_MONTHS) %>%
  distinct(PATIENT_ID, .keep_all = TRUE)

dim(clinical_unique) # Ska motsvara antal unika patienter


# ── 2. KONVERTERA OS_STATUS & DFS_STATUS TILL NUMERISK (0/1) ─────────────────
# TCGA-format brukar vara "0:LIVING" / "1:DECEASED" — checka dina värden:
table(clinical_unique$OS_STATUS)
#table(clinical_unique$DFS_STATUS)

# Konverterar till binär: 1 = event (deceased/recurrence), 0 = censurerad
clinical_unique$OS_event  <- ifelse(grepl("1|DECEASED|Dead",  clinical_unique$OS_STATUS,  ignore.case = TRUE), 1, 0)
#clinical_unique$DFS_event <- ifelse(grepl("1|Recurred|Progressed", clinical_unique$DFS_STATUS, ignore.case = TRUE), 1, 0)


# ── 3. KAPLAN-MEIER: OVERALL SURVIVAL BY SEX ─────────────────────────────────

fit_os_sex <- survfit(
  Surv(OS_MONTHS, OS_event) ~ SEX,
  data = clinical_unique
)

km_os_sex <- ggsurvplot(
  fit_os_sex,
  data          = clinical_unique,
  pval          = TRUE,          # Visar log-rank p-värde
  pval.method   = TRUE,          # Visar vilken test som används
  conf.int      = TRUE,          # 95% konfidensintervall
  risk.table    = TRUE,          # Antal i risk under grafen
  risk.table.height = 0.25,
  palette       = c("#E64B35", "#4DBBD5"),  # Röd = Female, Blå = Male
  legend.labs   = c("Female", "Male"),
  legend.title  = "Sex",
  xlab          = "Time (months)",
  ylab          = "Overall survival probability",
  title         = "Overall Survival by Sex in GBM",
  ggtheme       = theme_classic()
)

km_os_sex


# ── 4. KAPLAN-MEIER: DISEASE-FREE SURVIVAL BY SEX ────────────────────────────

fit_dfs_sex <- survfit(
  Surv(DFS_MONTHS, DFS_event) ~ SEX,
  data = clinical_unique
)

km_dfs_sex <- ggsurvplot(
  fit_dfs_sex,
  data          = clinical_unique,
  pval          = TRUE,
  pval.method   = TRUE,
  conf.int      = TRUE,
  risk.table    = TRUE,
  risk.table.height = 0.25,
  palette       = c("#E64B35", "#4DBBD5"),
  legend.labs   = c("Female", "Male"),
  legend.title  = "Sex",
  xlab          = "Time (months)",
  ylab          = "Disease-free survival probability",
  title         = "Disease-Free Survival by Sex in GBM",
  ggtheme       = theme_classic()
)

km_dfs_sex


# ── 5. KAPLAN-MEIER: OS STRATIFIERAT PÅ IMMUNCELLABUDANS (VALFRI CELLTYP) ────
# Byt ut "B_cell" mot valfri celltyp i ditt dataset, t.ex. "CD8_T_cell"

unique(data_immunoabundance$CELL_TYPE)
target_cell <- "Mast_cells"  # <-- Ändra till önskad celltyp

# Hämta abundansen för vald celltyp per patient
cell_abund <- data_immunoabundance %>%
  filter(CELL_TYPE == target_cell) %>%
  select(PATIENT_ID, ABUNDANCE)

# Merga med klinisk data
clinical_cell <- merge(clinical_unique, cell_abund, by = "PATIENT_ID")

# ── Bestäm optimal cutoff baserat på survival ───────────────────────────────
cutpoint <- surv_cutpoint(
  clinical_cell,
  time = "OS_MONTHS",
  event = "OS_event",
  variables = "ABUNDANCE"
)

clinical_cell_grouped <- surv_categorize(cutpoint)

table(clinical_cell_grouped$ABUNDANCE) # kontrollera gruppstorlek
cutpoint$cutpoint                      # visar vilken cutoff som valdes

fit_os_cell <- survfit(
  Surv(OS_MONTHS, OS_event) ~ ABUNDANCE,
  data = clinical_cell_grouped
)

km_os_cell <- ggsurvplot(
  fit_os_cell,
  data          = clinical_cell_grouped,
  pval          = TRUE,
  pval.method   = TRUE,
  conf.int      = TRUE,
  risk.table    = TRUE,
  risk.table.height = 0.25,
  palette       = c("#F39B7F", "#3C5488"),  # Orange = High, Mörkblå = Low
  legend.labs   = c("High abundance", "Low abundance"),
  legend.title  = target_cell,
  xlab          = "Time (months)",
  ylab          = "Overall survival probability",
  title         = paste("Overall Survival by", target_cell, "Abundance in GBM"),
  ggtheme       = theme_classic()
)

km_os_cell

# ── 7. LOG-RANK TEST SUMMARIES ────────────────────────────────────────────────

# Overall survival by sex
survdiff(Surv(OS_MONTHS, OS_event) ~ SEX, data = clinical_unique)

# Disease-free survival by sex
#survdiff(Surv(DFS_MONTHS, DFS_event) ~ SEX, data = clinical_unique)

# OS by immune cell abundance
survdiff(Surv(OS_MONTHS, OS_event) ~ ABUNDANCE, data = clinical_cell_grouped)
