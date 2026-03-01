library(biomaRt)  

rna <- read.delim("mrnadata_with_gene_names.txt", header=TRUE, stringsAsFactors=FALSE)

listEnsembl(
  mart = NULL,
  version = NULL,
  GRCh = NULL,
  mirror = NULL,
  verbose = FALSE
)

mart = useEnsembl('genes')
listDatasets(mart)

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
genes <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name","transcript_biotype"), filters = c("transcript_biotype","chromosome_name"),values = list("protein_coding",c(1:22)), mart = mart)


shouldImport=TRUE
saveFile="proteinCodingMouseGenes.Rda"
if (!shouldImport || !file.exists(saveFile)){
  print("Querying Biomart for protein coding genes")
  ensembl=useMart("ensembl")
  ensemblMouse = useDataset("mmusculus_gene_ensembl",mart=ensembl)
  mouseProteinCodingGenes = getBM(attributes=c("ensembl_gene_id","external_gene_name","description"), filters='biotype', values=c('protein_coding'), mart=ensemblMouse)
  save(mouseProteinCodingGenes,file=saveFile)
} else {
  print("Loading genes from savefile")
  load(saveFile)
}