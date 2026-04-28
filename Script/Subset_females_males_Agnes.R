getwd()
# Ladda ner clinical data och kolla på den
clin <- read.delim("AgnesVärld/Data/gbm_tcga_gdc/data_clinical_patient.txt",
                   header = TRUE,
                   comment.char = "#",
                   stringsAsFactors = FALSE,
                   check.names = FALSE)

clin <- clin[, c("PATIENT_ID", "SEX")]
head(clin)
View(clin)

# Ladda ner mrna data och kolla på den
mrna <- read.delim("Data/data_mrna_seq_tpm_proteinCoding.txt",
                   header = TRUE,
                   stringsAsFactors = FALSE,
                   check.names = FALSE)
View(mrna)

# Tar bort första kolumnen som innehåller gensymboler (kolumn 1)
patient_id_long <- colnames(mrna)[-1]
length(patient_id_long) # en kortare än mrna 
head(patient_id_long)

# Konvertera sample ID till pariten ID genom att ta bort de 4 sista karaktärerna. Totalt blir det 12
patient_id <- substr(patient_id_long, 1, 12)

colnames(mrna)[-1] <- patient_id
View(mrna)

# Skapar en nyckel med sample ID och patient ID för att kontrollera att allt stämmer
expr_key <- data.frame(
  PATIEN_ID_L  = patient_id_long,
  PATIENT_ID = patient_id,
  stringsAsFactors = FALSE
)
head(expr_key)

# Skapar en ny nyckel som mergar den gamla nyckeln med clinical data som filtrerats för sex och patient_ID
expr_key2 <- merge(expr_key, clin, by = "PATIENT_ID", all.x = TRUE) # all.x = TRUE innebär att all mrna-data sparas även om sex saknas
head(expr_key2)
table(expr_key2$SEX, useNA = "ifany") # Notera att det inte finns några NA!!!

head(expr_key2)

# Skapar en look_up-vektror där namnet är ID och värdet är könet eller sex 
sex_lookup <- setNames(expr_key2$SEX, expr_key2$PATIENT_ID)
sex_lookup["TCGA-02-0003"] # Exempel på hur det funkar 

# Ser till att orningen på alla kolumner matchar med den i nyckel 2 för att sedan kunna koppla
mrna2 <- mrna[, c("Gene_Symbol", expr_key2$PATIENT_ID)]
mrna2 <- colnames(mrna2)[-1]
head(mrna2)
View(mrna2)

# Omformaterar nyckeln till metadata och se till att Pateint ID matchar ordningen 
sample_meta <- expr_key2[, c("PATIENT_ID", "SEX")]
all(colnames(mrna2)[-1] == sample_meta$PATIENT_ID) # allt matchar!!!
View(sample_meta)

# Delar upp males och females i metadatan kopplat till ID
female_samples <- sample_meta$PATIENT_ID[sample_meta$SEX == "Female"]
male_samples   <- sample_meta$PATIENT_ID[sample_meta$SEX == "Male"]

length(female_samples) # Detta matchar med tabellen innan!
length(male_samples)

# Subsettar mrna tabellen med females rsp. males för att skapa två nya data-set 
mrna_female <- mrna[, c("Gene_Symbol", female_samples)]
mrna_male <- mrna[, c("Gene_Symbol", male_samples)]

dim(mrna_female) #dimensioner matchar med length då den ska vara 1 + antal kolumner 
dim(mrna_male)

# Skriver över tabellerna till nya txt-filer för att dela med gruppen
getwd()
write.table(mrna_female, #female mrna expression
            file = "Data/data_mrna_seq_tpm_female.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

write.table(mrna_male, #male gene expression 
            file = "Data/data_mrna_seq_tpm_male.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)