library(dplyr)

#### Candidate and related gene lists ####
# TODO : Change here if GENE LIST update
candidate_genes_path <-"../GENE_LIST/GENELIST_0913.txt"
candidate_genes <- read.table(candidate_genes_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

top_candidate_genes <- c("DSCAM", "KCNJ6", "RCAN1", "COL6A1", "COL6A2", "ALK2", "CRELD1", "CRELD2", "HEY2")
top_candidate_related_genes <- read.table("../GENE_LIST/target_gene_related.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)


#### GO and KEGG pathway gene lists #### 
cilium <- read.table("../GENE_LIST/GO0005929_cilium.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
protein_folding <- read.table("../GENE_LIST/GO0140662_ATPdependent protein folding chaperone.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# KEGG
ECM_interaction <- read.table("../GENE_LIST/hsa04512_ECM_receptor_interaction.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
hsa04310_Wnt <- read.table("../GENE_LIST/hsa04310_Wnt.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
hsa04330_Notch <- read.table("../GENE_LIST/hsa04330_Notch.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
hsa04340_hedgedog <- read.table("../GENE_LIST/hsa04340_hedgedog.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
HSA189451_heme <- read.table("../GENE_LIST/HSA189451_heme.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# STRING
CL8946_folic_acid <- read.table("../GENE_LIST/CL8946_folic_acid.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

#### HPO gene list ####
HP0006695 <- read.table("../GENE_LIST/genes_for_HP_0006695", header = TRUE, sep = "\t", stringsAsFactors = FALSE)$name
HP0006695 <- trimws(HP0006695)

#### To see if genes in genelist ####
# FOR TXT FILE
gene_list <- function(input, genelist) {
  # Ensure 'Gene_refgene' column exists (for Mviewer output, MAF<0.01)
  if (!"Gene_refgene" %in% colnames(input)) {
    stop("The input data does not contain a column named 'Gene_refgene'.")
  }
  
  # Filter based on the genelist
  output <- input %>% 
    filter(Gene_refgene %in% genelist$V1)
  
  return(output)
}


# N1675
top_candidate_genes_N1675_MAF0.01 = N1675_MAF0.01 %>% filter(Gene_refgene %in% top_candidate_genes)
top_candidate_related_genes_N1675_MAF0.01 = gene_list(N1675_MAF0.01,top_candidate_related_genes)
candidate_genes_N1675_MAF0.01 = gene_list(N1675_MAF0.01,candidate_genes)
ECM_interaction_N1675_MAF0.01 = gene_list(N1675_MAF0.01,ECM_interaction)
hsa04310_Wnt_N1675_MAF0.01 = gene_list(N1675_MAF0.01,hsa04310_Wnt)
cilium_N1675_MAF0.01 = gene_list(N1675_MAF0.01,cilium)
hsa04330_Notch_N1675_MAF0.01 = gene_list(N1675_MAF0.01,hsa04330_Notch)
hsa04310_Wnt_N1675_MAF0.01 = gene_list(N1675_MAF0.01,hsa04310_Wnt)
hsa04340_hedgedog_N1675_MAF0.01 = gene_list(N1675_MAF0.01,hsa04340_hedgedog)
CL8946_folic_acid_N1675_MAF0.01 = gene_list(N1675_MAF0.01,CL8946_folic_acid)
HSA189451_heme_N1675_MAF0.01 = gene_list(N1675_MAF0.01,HSA189451_heme)
HP0006695_N1675_MAF0.01 = N1675_MAF0.01 %>% filter(Gene_refgene %in% HP0006695)

# AI3008
top_candidate_genes_AI3008_MAF0.01 = AI3008_MAF0.01 %>% filter(Gene_refgene %in% top_candidate_genes)
top_candidate_related_genes_AI3008_MAF0.01 = gene_list(AI3008_MAF0.01,top_candidate_related_genes)
candidate_genes_AI3008_MAF0.01 = gene_list(AI3008_MAF0.01,candidate_genes)
ECM_interaction_AI3008_MAF0.01 = gene_list(AI3008_MAF0.01,ECM_interaction)
hsa04310_Wnt_AI3008_MAF0.01 = gene_list(AI3008_MAF0.01,hsa04310_Wnt)
cilium_AI3008_MAF0.01 = gene_list(AI3008_MAF0.01,cilium)
hsa04330_Notch_AI3008_MAF0.01 = gene_list(AI3008_MAF0.01,hsa04330_Notch)
hsa04310_Wnt_AI3008_MAF0.01 = gene_list(AI3008_MAF0.01,hsa04310_Wnt)
hsa04340_hedgedog_AI3008_MAF0.01 = gene_list(AI3008_MAF0.01,hsa04340_hedgedog)
CL8946_folic_acid_AI3008_MAF0.01 = gene_list(AI3008_MAF0.01,CL8946_folic_acid)
HSA189451_heme_AI3008_MAF0.01 = gene_list(AI3008_MAF0.01,HSA189451_heme)
HP0006695_AI3008_MAF0.01 = AI3008_MAF0.01 %>% filter(Gene_refgene %in% HP0006695)

