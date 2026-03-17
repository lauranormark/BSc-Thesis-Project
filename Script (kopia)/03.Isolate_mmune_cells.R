getwd()

#install.packages("remotes")
library(remotes)
#install.packages("devtools")
#devtools::install_github("cansysbio/ConsensusTME")
library(ConsensusTME)

#Laddar ner mrna-data med proteinkodande gener 
mrna_proteins <- read.delim(
  "Data/data_mrna_seq_tpm_proteinCoding.txt",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

#konvertera om till matris och se till att rad-namnen är gensymbolerna
expr_mat <- as.matrix(mrna_proteins[, -1]) # tar bort första kolumen (gensymboler)
rownames(expr_mat) <- mrna_proteins$Gene_Symbol # döper om radnamn i matrisen 
mode(expr_mat) <- "numeric"

#checkar för duplikat och na igen 
sum(duplicated(rownames(expr_mat)))
sum(is.na(rownames(expr_mat)))

dim(expr_mat) #dimension stämmer

#Konverterar om till log på expressionsdata eftersom paketet söker det 
expr_log2 <- log2(expr_mat + 1)

# Använder paketet för cancer GBM. Skapar matris med sammanställda värden från 18 immunceller samt immunscore
immune_res <- ConsensusTME::consensusTMEAnalysis(
  expr_log2,
  cancer = "GBM"
)
dim(immune_res) #Stämmer med 295 patienter och 19 rader totalt

write.table(
  immune_res,
  file = "Data/mrna_immune_cell_abundance.txt",
  sep = "\t",
  quote = FALSE
)
