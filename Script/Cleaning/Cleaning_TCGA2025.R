# Import dataset
setwd('Data') #ställ in wd för där din data finns, i dett fall var det script
tcga2025 <- read.table("gbm_tcga_pub2025_clinical_data.tsv", header=TRUE, sep="\t")

View(tcga2025)
summary(tcga2025)
table(tcga2025$Sex) # Undersök vilka variabler som finns per kolumn

# Drop unnecessairy variable, se nedan
drop <- c("Study.ID",
          "Patient.ID",
          "Sample:ID",
          "Biopsy.Site",
          "Cancer.Type",
          "Cancer.Type.Detailed",
          "Last.Communication.Contact.from.Initial.Pathologic.Diagnosis.Date",
          "Disease.Free..Months.",
          "Disease.Free.Status",
          "Ethnicity.Category",
          "Oncotree.Code",
          "Overall.Survival.Status",
          "Other.Patient.ID",
          "Other.Sample.ID",
          "Primary.Diagnosis",
          "Patient.Primary.Tumor.Site",
          "Prior.Malignancy",
          "Prior.Treatment",
          "Project.Identifier",
          "Project.Name",
          "Project.State",
          "Race.Category",
          "Number.of.Samples.Per.Patient",
          "Sample.Type",
          "Sample.type.id"
)

tcga2025_cleaned = tcga2025[,!(names(tcga2025) %in% drop)]
View(NAMNDATA_cleaned) # Noter att rengjord data är "NAMNDATA_cleaned"

# Check sex correspondance to variable. I detta fall CIMP
table(NAMNDATA_cleaned$Sex,
      NAMNDATA:cleaned$G.CIMP..Methylation)

# Split data for males and females
NAMNDATA_male <- NAMNDATA_cleaned[NAMNDATA_cleaned$Sex == "Male",]
NAMNDATA_female <- NAMNDATA_cleaned[NAMNDATA_cleaned$Sex == "Female",]

View(NAMNDATA_female)
View(NAMNDATA_male)