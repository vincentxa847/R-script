library(dplyr)
library(ggplot2)
library(ggrepel)

# -----------------------------------------------------------------------------
# Perform rare variant collapsing analysis
# GeneList relevant to phenotype and qualifying variants table of case and control will be input
# Create gene-by-individual collapsing matrix which indicates the presence of at least one QV
# Each gene is tested for an association between QV status and the phenotype of interest using Fisher Exact Test
# Fisher's Exact Test provides an exact p-value based on the observed data, regardless of sample size
# This is particularly beneficial when dealing with small datasets where the distribution of the data does not approximate normality
# Evaluate result by quantile-quantile plot (QQ plot)
# Reference : Rare- variant collapsing analyses for complex traits: guidelines and applications (Nat Rev Genet. 2019)
# -----------------------------------------------------------------------------
perform_fisher_test <- function(genelist, case_tables, control_tables) {
  
  # Helper function to create gene by individual collapsing matrix for cases or controls
  gene_by_individual_collapsing_matrix <- function(genelist, CaseorControl) {
    
    # Initialize an empty matrix with rows as genelist and columns as number of samples (tables in CaseorControl)
    result_matrix <- matrix(0, nrow = length(genelist), ncol = length(CaseorControl),
                            dimnames = list(genelist, paste0("Sample_", seq_along(CaseorControl))))
    
    # Loop over each sample (table) in CaseorControl
    for (i in seq_along(CaseorControl)) {
      # Get the current sample's table
      sample_table <- CaseorControl[[i]]
      
      # Filter the table by CADD_phred > 15 and Gene.refgene in the genelist
      filtered_table <- subset(sample_table, CADD_phred != "." & 
                                 round(as.numeric(CADD_phred), 2) >= 15  & Gene.refgene %in% genelist)
      
      # Count the number of occurrences for each gene in genelist
      gene_counts <- table(filtered_table$Gene.refgene)
      
      # Match gene names from gene_counts to row names of the result_matrix
      matching_genes <- match(names(gene_counts), rownames(result_matrix))
      
      # Update the result matrix with 1 (presence of QV) if gene is present, otherwise leave as 0
      result_matrix[matching_genes[!is.na(matching_genes)], i] <- 1
    }
    
    return(result_matrix)
  }
  
  # Create gene-by-individual matrices for cases and controls
  case_matrix <- gene_by_individual_collapsing_matrix(genelist, case_tables)
  control_matrix <- gene_by_individual_collapsing_matrix(genelist, control_tables)
  
  # Collapse rare variants per gene (sum across all individuals in cases and controls)
  case_collapsed <- rowSums(case_matrix)
  control_collapsed <- rowSums(control_matrix)
  
  # Number of cases and controls
  num_cases <- length(case_tables)
  num_controls <- length(control_tables)
  
  # Fisher's Exact Test for each gene
  p_values <- numeric(length(genelist))
  
  for (i in seq_along(genelist)) {
    # Create a 2x2 table for the current gene
    contingency_table <- matrix(c(
      case_collapsed[i],                   # Rare variants in cases
      num_cases - case_collapsed[i],      # No rare variant in cases
      control_collapsed[i],                # Rare variants in controls
      num_controls - control_collapsed[i]  # No rare variant in controls
    ), nrow = 2, byrow = TRUE)
    
    # Perform Fisher's Exact Test
    test <- fisher.test(contingency_table)
    p_values[i] <- test$p.value
  }
  
  # Return the p-values
  return(p_values)
}

# Example usage
genelist <- gene_lists$cilium
case_tables <- list(D25029_GeneList$cilium, D25046_GeneList$cilium)
control_tables <- list(D25007_GeneList$cilium)

# Perform the Fisher test with the original setup
p_values <- perform_fisher_test(genelist, case_tables, control_tables)

### Create a QQ-Plot
qqplot_pvalues <- function(pvals, genelist) {
  n <- length(pvals)
  
  # Expected values
  expected <- -log10((1:n) / n)
  
  # Observed values
  observed <- -log10(sort(pvals))
  
  # Get the indices of the top 5 smallest p-value
  top5_indices <- order(pvals)[1:5]
  
  # Get the top 5 gene names
  top5_genes <- genelist[top5_indices]
  
  # Calculate genomic inflation factor (lambda)
  lambda <- round(median(qchisq(pvals, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1), 3)
  
  # Calculate the 95% confidence intervals for the expected values (2.5th and 97.5th percentiles)
  ci_lower <- -log10(qbeta(0.025, seq(1, n), rev(seq(n, 1))))
  ci_upper <- -log10(qbeta(0.975, seq(1, n), rev(seq(n, 1))))
  
  # Prepare the data frame for ggplot2
  qq_data <- data.frame(
    expected = expected,
    observed = observed,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    is_top5 = FALSE  # Placeholder for top 5
  )
  
  # Mark the top 5 genes in the data
  qq_data$is_top5[top5_indices] <- TRUE
  
  # Create the QQ plot with ggplot2
  plot <- ggplot(qq_data, aes(x = expected, y = observed)) +
    geom_point(size = 1.5, color = "#1E4C9C") +  # Add the points for observed vs expected
    geom_line(aes(y = ci_lower), color = "#FBDFE2") +  # 2.5th percentile
    geom_line(aes(y = ci_upper), color = "#B83945") +  # 97.5th percentile
    labs(
      title = paste0("QQ-Plot Lambda = ", lambda),
      x = "Expected -log10(p-value)",
      y = "Observed -log10(p-value)"
    ) +
    theme_minimal() +  # Clean minimal theme
    theme(
      plot.title = element_text(hjust = 0.5),  # Center the plot title
      panel.grid = element_blank(),
      legend.position = "top"
    )
  
  return(plot)
}

# Plot the QQ-plot
qqplot_pvalues(p_values, gene_lists$cilium)
