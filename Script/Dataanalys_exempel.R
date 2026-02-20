# Import dataset
Glioblastoma_2013 <- read.table("gbm_tcga_pub2013_clinical_data.tsv", header=TRUE, sep="\t")
View(Glioblastoma_2013)

# Drop unnecessairy variables
drop <- c("Study.ID","Patient.ID", "Sample.ID", "Number.of.Samples.Per.Patient", "Sample.Type", "Somatic.Status", "therapy")
Glioblastoma_2013_cleaned = Glioblastoma_2013[,!(names(Glioblastoma_2013) %in% drop)]
View(Glioblastoma_2013_cleaned)

# Check sex correspondance to CIMP
table(Glioblastoma_2013$Sex,
      Glioblastoma_2013$G.CIMP..Methylation)

# split data for males and females
data_male <- Glioblastoma_2013_cleaned[Glioblastoma_2013_cleaned$Sex == "Male",]
data_female <- Glioblastoma_2013_cleaned[Glioblastoma_2013_cleaned$Sex == "Female",]

View(data_female)
View(data_male)
