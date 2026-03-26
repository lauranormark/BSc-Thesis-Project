getwd()
#Installera paket som är nödvändig
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("tidyr")
install.packages("dplyr")

#Läs in datan som data frame 
immune <- read.delim(
  "Filis/Datafiler/mrna_immune_cell_abundance.txt",
  header = TRUE,
  check.names = FALSE
)

dim(immune) #Samma som immunte_res

#Tar bort de 4 sista symbolerna i det långa patient-id för att matcha id med den kliniska datan
colnames(immune) <- substr(colnames(immune), 1, 12)

# Checkar för att se om det finns duplikat i patient ID. Dessa försvinner vid melt. 
truncated_ids <- substr(colnames(immune), 1, 12)
sum(duplicated(truncated_ids))
truncated_ids[duplicated(truncated_ids)] # visar exakt vilka som är duplikerade (9 st). 295 pat -> 286 pat

#konverar så att radnamn är nu en kolumn för att sedam kunna konverta om till lång data 
immune$CELL_TYPE <- rownames(immune)
rownames(immune) <- NULL

immune_long <- reshape2::melt(
  immune,
  id.vars = "CELL_TYPE",
  variable.name = "PATIENT_ID", 
  value.name = "ABUNDANCE"
  )

head(immune_long) # checkar så det ser rätt ut 
dim(immune_long)

#lägg till clinical data för att kunna jämföra
clinical <- read.delim("Filis/data_clinical_patient.txt",
                       header = TRUE, 
                       comment.char = "#", 
                       stringsAsFactors = FALSE, 
                       check.names = FALSE)

colnames(clinical)

#Loadar clinical data med värden som vi ska kolla på. 
clinical <- clinical[, c("PATIENT_ID",
                         "AGE",
                         "SEX",
                         "OS_STATUS",
                         "OS_MONTHS",
                         "DISEASE_TYPE",
                         "DAYS_TO_DEATH",
                         "PRIOR_TREATMENT")]


# Mergar med immune_long
immune_long <- merge(
  immune_long,
  clinical,
  by = "PATIENT_ID"
)
colnames(immune_long) # allt ser bra ut 

#checkar så att allt ser rätt ut 
dim(immune_long)
table(immune_long$SEX)
table(immune_long$CELL_TYPE) # totatl 286 patienter

#ser till att abundance är numerisk igen
immune_long$ABUNDANCE <- as.numeric(immune_long$ABUNDANCE)

write.table(immune_long, #male gene expression 
            file = "Filis/Datafiler/merged_immunoabundance_clinical.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

