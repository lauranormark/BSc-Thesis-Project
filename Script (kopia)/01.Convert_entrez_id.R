getwd() #Ska vara i ens egna solodev 
mrna <- read.delim("data_mrna_seq_tpm.txt",
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

mrna_geneSymbol <- mrna

# Konverterar entrez ID till gensymboler mha Geneconductor
mrna_geneSymbol$Gene_Symbol <- mapIds(org.Hs.eg.db,
                          keys = mrna$Entrez_Gene,
                          keytype = "ENTREZID",
                          column = "SYMBOL",
                          multiVals = "first")
# Nu har mrna data-setet fått en ny kolumn som heter Gene_Symbol med gensymboler

# Checka resultatet på nya data-setet
head(mrna_geneSymbol[, c("Entrez_Gene", "Gene_Symbol")])
sum(is.na(mrna_geneSymbol$Gene_Symbol))
sum(duplicated(mrna_geneSymbol$Gene_Symbol))

#ser att det finns 30 NA och 63 duplikat - tar bort dem 
mrna_geneSymbol_clean <- mrna_geneSymbol[!is.na(mrna_geneSymbol$Gene_Symbol) & !duplicated(mrna_geneSymbol$Gene_Symbol), ]

# Ta bort entrez gene ID kolumner 
mrna_geneSymbol_clean$Entrez_Gene <- NULL
mrna_geneSymbol_clean$Entrez_Gene_Id <- NULL

# Flyttar Gene_symbol först 
mrna_geneSymbol_clean <- mrna_geneSymbol_clean[, c("Gene_Symbol", setdiff(names(mrna_geneSymbol_clean), "Gene_Symbol"))]
which(names(mrna_geneSymbol_clean) == "Gene_Symbol") # checkar att allt gick bra

# Laddar ner datan igen i ett nytt data-set som heter samma sak fast med _geneSymbol 
write.table(mrna_geneSymbol_clean,
            file = "Datafiler/data_mrna_seq_tpm_geneSymbols.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
