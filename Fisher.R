library(dplyr)
#### Prepare data for Fisher test ####
## Function to count number of variant in each gene of gene list
variants_in_genes <- function(target_gene, variants, name){
  # count variant number of each genes and set column name
  raw <- variants %>%
    group_by(Gene_refgene) %>%
    summarise(!!sym(name) := n())  
  
  # extract data of gene in gene list and replace NA with 0 for missing counts
  gene_and_count <- target_gene %>%
    left_join(raw, by = "Gene_refgene") %>%
    mutate(!!sym(name) := ifelse(is.na(!!sym(name)), 0, !!sym(name)))
  
  return (gene_and_count)
}

## Function to apply the function to each sample in each group and store the results
combine_variants_in_genes <- function(group, gene_list, sample_names) {
  results_list <- list()
  
  # Loop through each group sample and apply the variants_in_genes function
  for (i in seq_along(group)) {
    variants <- group[[i]]
    name <- sample_names[i]
    
    results_list[[i]] <- variants_in_genes(gene_list, variants, name)
  }
  
  # Combine the results into one data frame using full_join
  final_results <- Reduce(function(x, y) full_join(x, y, by = "Gene_refgene"), results_list)
  
  # View the final combined result
  print(final_results)
  
  # Return the combined result
  return(final_results)
}

#### Fisher’s exact test on mutation burden ####
Fisher_mutation_burden <- function(case, control) {
  
  # Count how many samples in the case group have at least one variant
  samples_with_variants_case <- sum(rowSums(as.data.frame(case[,-1])) > 0)  # Count of samples with at least one variant
  total_case_samples <- nrow(case)  # Total number of samples in case group
  
  # Count how many samples in the control group have at least one variant
  samples_with_variants_control <- sum(rowSums(as.data.frame(control[,-1])) > 0)  # Count of samples with at least one variant
  total_control_samples <- nrow(control)  # Total number of samples in control group
  
  # Calculate the number of samples without variants
  samples_without_variants_case <- total_case_samples - samples_with_variants_case
  samples_without_variants_control <- total_control_samples - samples_with_variants_control
  
  # Create a 2x2 contingency table
  table <- matrix(c(
    samples_with_variants_case, samples_without_variants_case,  # case group: presence vs absence of variants
    samples_with_variants_control, samples_without_variants_control  # control group: presence vs absence of variants
  ), nrow = 2)
  
  # Perform Fisher's exact test
  test <- fisher.test(table)
  
  # Return the p-value, odds ratio, and confidence intervals
  results <- data.frame(
    p_value = test$p.value,
    OR_fisher = test$estimate,
    OR_fisher_ci_left = test$conf.int[1],
    OR_fisher_ci_right = test$conf.int[2]
  )
  
  return(results)
}

#### Analysis of the mutation burden of rare non-synonymous variants (missense) in cilium components across DS-ECD and non-DS-ECD samples ####
## Gene List to analysis (cilium_components in this case)
cilium_components <- read.csv("../GENE_LIST/cilium_components.txt",
                              header = FALSE, sep = "\t", 
                              stringsAsFactors = FALSE, col.names = "Gene_refgene")


## Combine rare non-synonymous variants of samples of each group into a list 
DS_ECD <- list(D25029_MAF0.01_missense, D25046_MAF0.01_missense)
nonDS_ECD <- list(D25007_MAF0.01_missense)
# TODO: Change here
sample_names <- c("D25007_MAF0.01_missense")

Fisher_cilium_components_DSECD <- combine_variants_in_genes(DS_ECD, cilium_components, sample_names)
Fisher_cilium_components_nonDSECD <- combine_variants_in_genes(nonDS_ECD, cilium_components, sample_names)

test <- Fisher_mutation_burden(Fisher_cilium_components_DSECD,Fisher_cilium_components_nonDSECD)
