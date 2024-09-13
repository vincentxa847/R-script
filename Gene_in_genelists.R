library(dplyr)

#### Candidate and related gene lists ####
# TODO : Change here if GENE LIST update
candidate_genes_path <-"D:/DS Project_Vincent/Data/GENE_LIST/GENELIST_0913.txt"

candidate_genes <- read.table(candidate_genes_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
top_candidate_genes <- c("DSCAM", "KCNJ6", "RCAN1", "COL6A1", "COL6A2", "ALK2", "CRELD1", "CRELD2", "HEY2")
top_candidate_related_genes <- read.table("D:/DS Project_Vincent/Data/GENE_LIST/target_gene_related.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)


#### GO and KEGG pathway gene lists #### 
# GO 
cilium <- read.table("D:/DS Project_Vincent/Data/GENE_LIST/GO0005929_cilium.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
protein_folding <- read.table("D:/DS Project_Vincent/Data/GENE_LIST/GO0140662_ATPdependent protein folding chaperone.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# KEGG
ECM_interaction <- read.table("D:/DS Project_Vincent/Data/GENE_LIST/hsa04512_ECM_receptor_interaction.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
hsa04310_Wnt <- read.table("D:/DS Project_Vincent/Data/GENE_LIST/hsa04310_Wnt.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
hsa04330_Notch <- read.table("D:/DS Project_Vincent/Data/GENE_LIST/hsa04330_Notch.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
hsa04340_hedgedog <- read.table("D:/DS Project_Vincent/Data/GENE_LIST/hsa04340_hedgedog.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
HSA189451_heme <- read.table("D:/DS Project_Vincent/Data/GENE_LIST/HSA189451_heme.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# STRING
CL8946_folic_acid <- read.table("D:/DS Project_Vincent/Data/GENE_LIST/CL8946_folic_acid.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

#### To see if genes in genelist ####
gene_list <- function(input, genelist) {
  # Ensure 'Gene_refgene' column exists
  if (!"Gene_refgene" %in% colnames(input)) {
    stop("The input data does not contain a column named 'Gene_refgene'.")
  }
  
  # Filter based on the genelist
  output <- input %>% 
    filter(Gene_refgene %in% genelist)
  
  return(output)
}

cilium_N1675 = gene_list(N1675_MAF0.01,ECM_interaction)

