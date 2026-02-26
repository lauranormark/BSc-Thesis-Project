#läs in textfilen
rna <- read.delim("data_mrna_seq_tpm.txt", header=TRUE, stringsAsFactors=FALSE)

#ladda bioconductuctor paketen
library(org.Hs.eg.db)
library(AnnotationDbi)

#Ändrar entrez  
rna$Entrez_Gene_Id <- as.character(rna$Entrez_Gene_Id)

#skapa en nu kolumn med genen name
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(rna$Entrez_Gene_Id),
  keytype = "ENTREZID",
  columns = "SYMBOL"
)

gene_map <- gene_map[!duplicated(gene_map$ENTREZID), ]

rna$GENENAME <- gene_map$SYMBOL[
  match(rna$Entrez_Gene_Id, gene_map$ENTREZID)
]

rna <- rna[, c("Entrez_Gene_Id",  #flytta genename till kolumn 2
               "GENENAME", 
               setdiff(colnames(rna), 
                       c("Entrez_Gene_Id", 
                         "GENENAME")))]

#ta bort entrez kolumn samt spara ny fil
rna$Entrez_Gene_Id <- NULL
write.table(rna, "mrnadata_with_gene_names.txt", sep="\t", row.names=FALSE, quote=FALSE)
