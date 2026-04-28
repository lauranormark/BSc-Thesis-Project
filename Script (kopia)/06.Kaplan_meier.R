getwd()
#setwd("AgnesVärld")

#Install packages if needed:
#install.packages("survival")
#install.packages("survminer")

library(survival)
library(survminer)
library(dplyr)

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────

data_immunoabundance <- read.delim(
  "Data/merged_immunoabundance_clinical.txt",
  header = TRUE,
  check.names = FALSE
)

data_immunoabundance_female <- data_immunoabundance[data_immunoabundance$SEX == "Female", ]
data_immunoabundance_male   <- data_immunoabundance[data_immunoabundance$SEX == "Male", ]

# Kaplan-Meier behöver en rad per patient (inte long-format).
# Vi tar ut unika patienter med klinisk data.
clinical_unique <- data_immunoabundance %>%
  select(PATIENT_ID, AGE, SEX, OS_STATUS, OS_MONTHS) %>%
  distinct(PATIENT_ID, .keep_all = TRUE)

# for FEMALES
clinical_unique_female <- data_immunoabundance_female %>%
  select(PATIENT_ID, AGE, SEX, OS_STATUS, OS_MONTHS) %>%
  distinct(PATIENT_ID, .keep_all = TRUE)

# for MALES
clinical_unique_male <- data_immunoabundance_male %>%
  select(PATIENT_ID, AGE, SEX, OS_STATUS, OS_MONTHS) %>%
  distinct(PATIENT_ID, .keep_all = TRUE)

dim(clinical_unique) # Ska motsvara antal unika patienter
dim(clinical_unique_female)
dim(clinical_unique_male)

# ── 2. KONVERTERA OS_STATUS & DFS_STATUS TILL NUMERISK (0/1) ─────────────────
# TCGA-format brukar vara "0:LIVING" / "1:DECEASED" — checka dina värden:
table(clinical_unique$OS_STATUS)
table(clinical_unique_female$OS_STATUS)
table(clinical_unique_male$OS_STATUS)


# Konverterar till binär: 1 = event (deceased/recurrence), 0 = censurerad
clinical_unique$OS_event  <- ifelse(grepl("1|DECEASED|Dead",  clinical_unique$OS_STATUS,  ignore.case = TRUE), 1, 0)
clinical_unique_female$OS_event  <- ifelse(grepl("1|DECEASED|Dead",  clinical_unique_female$OS_STATUS,  ignore.case = TRUE), 1, 0)
clinical_unique_male$OS_event  <- ifelse(grepl("1|DECEASED|Dead",  clinical_unique_male$OS_STATUS,  ignore.case = TRUE), 1, 0)

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


# ── 5. KAPLAN-MEIER: OS STRATIFIERAT PÅ IMMUNCELLABUDANS (VALFRI CELLTYP), FEMALES ────
# Byt ut "B_cell" mot valfri celltyp i ditt dataset, t.ex. "CD8_T_cell"
target_cell <- "Immune_Score"  # <-- Ändra till önskad celltyp

# FEMALES
unique(data_immunoabundance_female$CELL_TYPE)
# Hämta abundansen för vald celltyp per patient
cell_abund_female <- data_immunoabundance_female %>%
  filter(CELL_TYPE == target_cell) %>%
  select(PATIENT_ID, ABUNDANCE)

# Merga med klinisk data
clinical_cell_female <- merge(clinical_unique_female, cell_abund_female, by = "PATIENT_ID")

# ── Bestäm optimal cutoff baserat på survival ───────────────────────────────
cutpoint_female <- surv_cutpoint(
  clinical_cell_female,
  time = "OS_MONTHS",
  event = "OS_event",
  variables = "ABUNDANCE"
)

clinical_cell_grouped_female <- surv_categorize(cutpoint_female)

table(clinical_cell_grouped_female$ABUNDANCE) # kontrollera gruppstorlek
cutpoint_female$cutpoint                      # visar vilken cutoff som valdes

fit_os_cell_female <- survfit(
  Surv(OS_MONTHS, OS_event) ~ ABUNDANCE,
  data = clinical_cell_grouped_female
)

km_os_cell_female <- ggsurvplot(
  fit_os_cell_female,
  data          = clinical_cell_grouped_female,
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
  title         = paste("Overall Survival by", target_cell, "Abundance in GBM (females)"),
  ggtheme       = theme_classic()
)

km_os_cell_female

#MALES
unique(data_immunoabundance_male$CELL_TYPE)

# Hämta abundansen för vald celltyp per patient
cell_abund_male <- data_immunoabundance_male %>%
  filter(CELL_TYPE == target_cell) %>%
  select(PATIENT_ID, ABUNDANCE)

# Merga med klinisk data
clinical_cell_male <- merge(clinical_unique_male, cell_abund_male, by = "PATIENT_ID")

# ── Bestäm optimal cutoff baserat på survival ───────────────────────────────
cutpoint_male <- surv_cutpoint(
  clinical_cell_male,
  time = "OS_MONTHS",
  event = "OS_event",
  variables = "ABUNDANCE"
)

clinical_cell_grouped_male <- surv_categorize(cutpoint_male)

table(clinical_cell_grouped$ABUNDANCE) # kontrollera gruppstorlek
cutpoint_male$cutpoint                     # visar vilken cutoff som valdes

fit_os_cell_male <- survfit(
  Surv(OS_MONTHS, OS_event) ~ ABUNDANCE,
  data = clinical_cell_grouped_male
)

km_os_cell_male <- ggsurvplot(
  fit_os_cell_male,
  data          = clinical_cell_grouped_male,
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
  title         = paste("Overall Survival by", target_cell, "Abundance in GBM (male)"),
  ggtheme       = theme_classic()
)

km_os_cell_male

# ── 7. LOG-RANK TEST SUMMARIES ────────────────────────────────────────────────

# Overall survival by sex
survdiff(Surv(OS_MONTHS, OS_event) ~ SEX, data = clinical_unique)

# OS by immune cell abundance
survdiff(Surv(OS_MONTHS, OS_event) ~ ABUNDANCE, data = clinical_cell_grouped)
