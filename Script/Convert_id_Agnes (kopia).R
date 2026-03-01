setwd("AgnesVärld") #Denna ska vara unik för varje användare 

mrna <- read.delim("Data/gbm_tcga_gdc/data_mrna_seq_tpm.txt",
                   header = TRUE,
                   stringsAsFactors = FALSE,
                   check.names = FALSE)

# Checka hur datasetet ser ut
str(mrna[1:5, 1:5]) # struktur och data-typ av de första 5 raderna i de första 5 kolumnerna 
head(mrna$Entrez_Gene) # checkar de första värdena i kolumnen Entrez_Gene

# Installera packeten som krävs för att konvertera entrez gene ID till gen-namn. Se till att BiocManager är installerat innan.
library(AnnotationDbi)
library(org.Hs.eg.db)

# Ser till att Entrez_Gene är karaktär
mrna$Entrez_Gene <- as.character(mrna$Entrez_Gene)

# Konverterar entrez ID till gensymboler mha Geneconductor
mrna$Gene_Symbol <- mapIds(org.Hs.eg.db,
                          keys = mrna$Entrez_Gene,
                          keytype = "ENTREZID",
                          column = "SYMBOL",
                          multiVals = "first")
# Nu har mrna data-setet fått en ny kolumn som heter Gene_Symbol med gensymboler

# Checka resultatet på nya data-setet
head(mrna[, c("Entrez_Gene", "Gene_Symbol")])
table(is.na(mrna$GeneSymbol))

# Ta bort entrez gene ID kolumner 
mrna$Entrez_Gene <- NULL
mrna$Entrez_Gene_Id <- NULL

# Flyttar Gene_symbol först 
mrna <- mrna[, c("Gene_Symbol", setdiff(names(mrna), "Gene_Symbol"))]
which(names(mrna) == "Gene_symbol") # checkar att allt gick bra

# Droppa rader med inga gen-symboler
mrna2 <- mrna[!is.na(mrna$Gene_Symbol), ]

# Laddar ner datan igen i ett nytt data-set som heter samma sak fast med _geneSymbol 
write.table(mrna2,
            file = "Data/data_mrna_seq_tpm_geneSymbols.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)