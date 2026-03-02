getwd()

# Installera bioaRt i denna fil och selekterar ensemble av gener som ska jömföras mot med mart
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Ladda ner datan med Gen-symboler (Inte Entrez gene ID)
mrna <- read.delim("Data/data_mrna_seq_tpm_geneSymbols.txt",
                   header = TRUE,
                   stringsAsFactors = FALSE,
                   check.names = FALSE)

head(mrna[1:5, 1:5]) #checkar att datan ser rätt ut

# Selektera kolumnen med alla gener: första kolumenen som heter Gene_Symbols
genes <- mrna$Gene_Symbol # selekterar första kolumnen
genes <- unique(genes) # tar bort duplikat av rader/element
genes <- genes[!is.na(genes) & genes != ""] # Tar bort NA och tomma rader/element

sum(duplicated(genes))

# Frågar Ensembl: “For this list of gene symbols, tell me their gene type (biotype).” och storar det i annot
annot <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype"), # Vad den frågar efter: gensymbol och dess biotyp
  filters    = "hgnc_symbol", # Vilken typ av gen-identifikation använder jag: symboler
  values     = genes, # Vilken vektor undersöker jag?
  mart       = mart # Specificerat ovan: filtrerar inom dataset = "hsapiens_gene_ensembl" eller mänskliga gener
)

head(annot) # Undersök hur annot ser ut

# Filtrerar för de gener som säger att de är proteinkodande samt ser till att de bara upprepas en gång
# Skapar en vektor med namnet på alla gener som är protein-kodande
protein_genes <- unique(annot$hgnc_symbol[annot$gene_biotype == "protein_coding"])
head(protein_genes) # visade sig vara totalt 19361 gener vilket stämmer med litteraturen
length(protein_genes)
View(protein_genes)

# Subsettar expressionsdatan med enbart de gener som kodar för proteiner 
mrna_protein <- mrna[mrna$Gene_Symbol %in% protein_genes, ]
nrow(mrna_protein)
ncol(mrna_protein)
View(mrna_protein)

# Undersök duplikat så att längden matchar för att ta bort duplikat 
length(protein_genes)
nrow(mrna_protein)
length(unique(mrna_protein$Gene_Symbol))
sum(duplicated(mrna_protein$Gene_Symbol))

# Ta bort duplikat 
mrna_protein_unique <- mrna_protein[!duplicated(mrna_protein$Gene_Symbol), ]
View(mrna_protein_unique)
nrow(mrna_protein_unique)

# Laddar ner datan igen i ett nytt data-set som heter samma sak fast med _proteinCoding 
getwd()
write.table(mrna_protein_unique,
            file = "Data/data_mrna_seq_tpm_proteinCoding.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)