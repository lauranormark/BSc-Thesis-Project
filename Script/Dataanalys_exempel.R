# Import dataset
Glioblastoma_2013 <- read.table("gbm_tcga_pub2013_clinical_data.tsv", header=TRUE, sep="\t")
View(Glioblastoma_2013)

# Drop unnecessairy variables
drop <- c("Study.ID","Patient.ID", "Sample.ID", "Number.of.Samples.Per.Patient", "Sample.Type", "Somatic.Status", "therapy")
Glioblastoma_2013_cleaned = Glioblastoma_2013[,!(names(Glioblastoma_2013) %in% drop)]
View(Glioblastoma_2013_cleaned)

# Group into males and females

