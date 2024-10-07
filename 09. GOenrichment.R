library("reshape2")
library("biomaRt")
library("topGO")

#### Fetch geneId2GO and Gene Universe of human gene ####
### geneId2GO
# Connect to the Ensembl database
ensembl <- useMart("ensembl")
# Set the dataset to human genes
ensembl_human <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

go_mapping <- getBM(attributes = c("ensembl_gene_id", "go_id"),
                    filters = "biotype",
                    values = "protein_coding",  # You can change this to other gene types if needed
                    mart = ensembl_human)

# View the first few rows of the GO mapping
head(go_mapping)

# Take some time to fetch, so store as rds
saveRDS(go_mapping, file = "go_mapping.rds")

### Gene Universe (all protein-coding genes) [TODO:須確定是否該用 all protein-coding genes, 目前找不到 GO 導致 ERROR: 收捲時發生錯誤: 第一引數必須是向量]
# Generate the gene universe as a unique list of Ensembl gene IDs
gene_universe <- unique(go_mapping$ensembl_gene_id)
# View the number of genes in the universe
length(gene_universe) # 23231


#### GOenrichment Analysis ####
GOenrichment = function(genesOfInterest) {

  # Create the geneList factor (1 = genes of interest, 0 = background genes)
  geneList <- factor(as.integer(gene_universe %in% gene_mapping$ensembl_gene_id))
  names(geneList) <- gene_universe
  
  head(geneList)
  
  ## Create the topGO object for Biological Process (BP) ontology
  BPenrichment <- new("topGOdata", 
                      description="Biological Process Enrichment Analysis", 
                      ontology="BP", 
                      allGenes=geneList, 
                      annot=annFUN.gene2GO, 
                      gene2GO=go_mapping, 
                      nodeSize=10)  # nodeSize set to filter out small GO terms
  
  ## Perform Fisher's test for enrichment
  BPFisherResults <- runTest(BPenrichment, statistic="fisher")
  
  ## Retrieve and filter enriched GO terms (only those with p-value <= 0.05)
  BPenriched <- GenTable(BPenrichment, 
                         Fisher=BPFisherResults, 
                         orderBy="Fisher", 
                         topNodes=20)  # Adjust topNodes if needed
  BPsignificant <- BPenriched[as.numeric(BPenriched$Fisher) <= 0.05, ]
  
  ## Return the table of significant GO terms
  return(BPsignificant)
}



#### Perform Analysis ####
genes_of_interest <- D25007_MAF0.01_CADD15$Gene_refgene

# Query Ensembl for Ensembl IDs based on genes_of_interest
gene_mapping <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                      filters = "external_gene_name",
                      values = genes_of_interest,
                      mart = ensembl_human)

# View the result
print(gene_mapping$ensembl_gene_id)
GO_ <- GOenrichment(gene_mapping$ensembl_gene_id)
