getwd()

# Installera bioaRt i denna fil och selekterar ensemble av gener som ska jömföras mot med mart
library(biomaRt)

mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

# Hej Fili//M.Cuz
# Ladda ner datan med Gen-symboler (Inte Entrez gene ID)
mrna_geneSymbol_clean <- read.delim("Datafiler/data_mrna_seq_tpm_geneSymbols.txt",
                   header = TRUE,
                   stringsAsFactors = FALSE,
                   check.names = FALSE)

head(mrna_geneSymbol_clean[1:5, 1:5]) #checkar att datan ser rätt ut

# Selektera kolumnen med alla gener: första kolumenen som heter Gene_Symbols
genes <- mrna_geneSymbol_clean$Gene_Symbol # selekterar första kolumnen
sum(is.na(genes)) # Vad är det här? :)
sum(duplicated(genes)) #Inga na eller duplikat

# Frågar Ensembl: “For this list of gene symbols, tell me their gene type (biotype).” och storar det i annot
library(biomaRt)

mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror  = "useast"
)

annot <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype"),
  filters    = "hgnc_symbol",
  values     = genes,
  mart       = mart
)

head(annot) # Undersök hur annot ser ut

# Filtrerar för de gener som säger att de är proteinkodande samt ser till att de bara upprepas en gång
# Skapar en vektor med namnet på alla gener som är protein-kodande
protein_genes <- unique(annot$hgnc_symbol[annot$gene_biotype == "protein_coding"])
length(protein_genes)  # visade sig vara totalt 19361 gener vilket stämmer med litteraturen

# Subsettar expressionsdatan med enbart de gener som kodar för proteiner 
mrna_proteins <- mrna_geneSymbol_clean[mrna_geneSymbol_clean$Gene_Symbol %in% protein_genes, ]
nrow(mrna_proteins)
ncol(mrna_proteins)

# Undersök duplikat så att längden matchar för att ta bort duplikat 
length(protein_genes)
nrow(mrna_proteins)

# Laddar ner datan igen i ett nytt data-set som heter samma sak fast med _proteinCoding 
write.table(mrna_proteins,
            file = "Datafiler/data_mrna_seq_tpm_proteinCoding.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
