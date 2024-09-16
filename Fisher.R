library(readr)
library(dplyr)
library(stats)

#### Data Preprocessing ####
# Data from MVIEWER (MAF<0.01)
N1675_rare <-  read.table("../N1675/N1675-MAF0.01-cleanup.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
AI3008_rare <-  read.table("../AI3008/AI3008-MAF0.01-cleanup.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Handle gene list
TOP_CANDITATE_RELATED_GENES <- read.table("../GENE_LIST/target_gene_related.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# make left join key same
colnames(TOP_CANDITATE_RELATED_GENES) <- c('Gene_refgene') 
TOP_CANDITATE_RELATED_GENES <- unique(TOP_CANDITATE_RELATED_GENES) # Remove duplicate rows

#### Count number of variant in each target gene in case and control group ####
# Count number of variant in each target gene
variants_in_genes <- function(target_gene, variants, name){
  # Step 1: Count occurrences of genes in case or control
  raw <- variants %>%
    group_by(Gene_refgene) %>%
    summarise(!!sym(name) := n())  # Dynamically set the column name
  
  # Step 2: Left join with target_gene to ensure all genes are included
  gene_and_count <- target_gene %>%
    left_join(raw, by = "Gene_refgene") %>%
    mutate(!!sym(name) := ifelse(is.na(!!sym(name)), 0, !!sym(name)))  # Replace NA with 0 for missing counts
  
  return (gene_and_count)
}

N1675_rare_counts <- variants_in_genes(TOP_CANDITATE_RELATED_GENES,N1675_rare,"N1675_rare")
AI3008_rare_counts <- variants_in_genes(TOP_CANDITATE_RELATED_GENES,AI3008_rare,"AI3008_rare")

# Group the tables of same group together
merge_variant_counts <- function(group_column, ..., count_columns) {
  tables <- list(...)  # Collect all the tables passed as arguments
  table_names <- count_columns  # Names for the count columns
  
  # Initialize the first table as the base for merging
  merged_table <- tables[[1]] %>%
    select(all_of(group_column), all_of(table_names[1]))  # Select group column and the first count column
  
  # Loop through the remaining tables and merge them
  for (i in 2:length(tables)) {
    merged_table <- merged_table %>%
      left_join(select(tables[[i]], all_of(group_column), all_of(table_names[i])), 
                by = group_column)  # Merge on the group column
  }
  
  # Set the group column as row names and remove it from the final table
  merged_table_with_row_names <- merged_table %>%
    tibble::column_to_rownames(var = group_column)
  
  return(merged_table_with_row_names)
}

# Merging example 
merged_table <- merge_variant_counts(
  group_column = "Gene_refgene",  # Column to merge on
  N1675_rare_counts, AI3008_rare_counts,  # Tables to merge
  count_columns = c("N1675_rare", "AI3008_rare")  # Column names for counts
)


#### Fisher’s exact test ####
Fisher <- function(target_gene,case,control) {
  # create dataframe for result
  results <- data.frame(gene = target_gene$Gene_refgene, p_value = NA, odd_ratio = NA) 
  
  # case and control should be in an iterable list
  # count the occurance of variants in each gene and add together 
  for (gene in target_gene$Gene_refgene) {
    case_count <- case[case$Gene_refgene == gene,"Number"]
    if (case_count != 0) { # if error : length > 1 means gene has dupication
      case_count = 1 # case has variant in this genes 
    } else {
      case_count = 0}
    
    control_count <- control[control$Gene_refgene == gene,"Number"]
    if (control_count != 0) {
      control_count  = 1 # case has variant in this genes 
    } else {
      control_count  = 0}
    
    # Create 2x2 Contingency Table for Fisher exact test
    table <- matrix(c(case_count, 1 - case_count,  # sum(case_filtered$Gene_refgene == gene) - case_count
                      control_count, 1 - control_count),
                    nrow = 2)
    test <- fisher.test(table)
    results$p_value[results$Gene_refgene == gene] <- test$p.value
    results$odd_ratio [results$Gene_refgene == gene] <- test$estimate  
  }
  return(results)
}

N1675_AI3008 <- Fisher(TOP_CANDITATE_RELATED_GENES,case_counts,count_counts)
