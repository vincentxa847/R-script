library(SKAT)

# -----------------------------------------------------------------------------
#### Prepare the Data for SKAT (Genotype Matrix, Phenotype Vector) #### 
# -----------------------------------------------------------------------------
## Function to extract genotype data for either a gene list or a single gene
extract_genotype <- function(sample_data, genes) {
  
  # Filter the sample data
  filtered_data <- sample_data[sample_data$Gene_refgene %in% genes, ]
  
  # Create a variant identifier by combining "Chr" and "Start"
  filtered_data$Variant <- paste(filtered_data$Chr, filtered_data$Start, sep = "_")
  
  # Convert genotypes ('het' -> 1, 'hom' -> 2)
  genotypes <- ifelse(filtered_data$Otherinfo1 == "het", 1, 2)
  
  # Return a data frame with variant identifier and genotypes
  return(data.frame(Variant = filtered_data$Variant, Genotype = genotypes))
}

## Function to create genotype matrix and phenotype vector for SKAT
generate_SKAT_data <- function(case_list, control_list, gene_list) {
  
  # Initialize lists to store genotype data and all observed variants
  genotype_list <- list()
  all_variants <- c() # Collect all unique variants
  sample_names <- c() # Collect all samples
  
  # Process case data
  for (i in seq_along(case_list)) {
    sample_genotype <- extract_genotype(case_list[[i]], gene_list)
    genotype_list[[paste0("Case_", i)]] <- sample_genotype
    all_variants <- union(all_variants, sample_genotype$Variant)  
    sample_names <- c(sample_names, paste0("Case_", i))
  }
  
  # Process control data
  for (i in seq_along(control_list)) {
    sample_genotype <- extract_genotype(control_list[[i]], gene_list)
    genotype_list[[paste0("Control_", i)]] <- sample_genotype
    all_variants <- union(all_variants, sample_genotype$Variant)  
    sample_names <- c(sample_names, paste0("Control_", i))
  }
  
  # If no variants are found, return NULL
  if (length(all_variants) == 0) return(NULL)
  
  # Create an empty genotype matrix with all variants as columns and samples as rows
  # default value is 0, which means homozygous for the major allele (no variant)
  genotype_matrix <- matrix(0, nrow = length(genotype_list), ncol = length(all_variants))
  colnames(genotype_matrix) <- all_variants
  rownames(genotype_matrix) <- sample_names
  
  # Fill in the genotype matrix with data
  for (i in seq_along(genotype_list)) {
    sample_genotype <- genotype_list[[i]]
    variant_indices <- match(sample_genotype$Variant, all_variants)  # Find column indices for this sample's variants
    genotype_matrix[i, variant_indices] <- sample_genotype$Genotype  # Fill in the genotypes
  }
  
  # Create phenotype vector (1 for cases, 0 for controls)
  phenotype_vector <- c(rep(1, length(case_list)), rep(0, length(control_list)))
  
  return(list(genotype_matrix = genotype_matrix, phenotype_vector = phenotype_vector))
}

## Main Function to Run SKAT by Gene List or by Single Gene
run_SKAT <- function(geneList, group1, group2, SKAT_or_SKATO, by_gene = FALSE){
  
  genelist2analysis <- geneList

  if (SKAT_or_SKATO == "SKAT") {
    kernel2use <- "davies"  # SKAT method
  } else if (SKAT_or_SKATO == "SKAT-O") {
    kernel2use <- "optimal.adj"  # SKAT-O method
  }
  
  # RUN SKAT BY GENELIST OR BY GENE
  if (by_gene) { 
    
    # list to hold result of each gene
    skat_results <-list()
    
    for (gene in genelist2analysis) {
      
      # Generate genotype matrix and phenotype vector for the gene
      skat_data <- generate_SKAT_data(group1, group2, gene)
      
      # Skip if no data was found for the gene
      if (is.null(skat_data)) {
        next
      }
      
      # Run SKAT
      # Create null model
      null_model <- SKAT_Null_Model(skat_data$phenotype_vector ~ 1, out_type = "D", Adjustment = TRUE)
      print(paste("Running SKAT for gene:", gene))
      set.seed(123)
      skat_result <- SKAT(skat_data$genotype_matrix, null_model, kernel = "linear.weighted", method = kernel2use, weights.beta = c(1, 25))
      
      # Store results for the gene
      skat_results[[gene]] <- skat_result
    }
    
  } else {
    
    ##  Define covariate (for now sample size is too small (n=3) to include covariate)
    # include sex as covariate (1 = Male, 0 = Female)
    #sex_vector <- c(1, 0, 0)  # D25029:M, D25046:F, D25007:F
    # Create a data frame for covariates
    #covariate_df <- data.frame(sex = sex_vector)
    
    ## Generate the genotype matrix and phenotype vector for SKAT analysis
    skat_data <- generate_SKAT_data(group1, group2, genelist2analysis)
    
    ## Run SKAT
    # Creating the Null Model, which represents the baseline against which the alternative hypothesis
    # Not using covariate here (~ 1), which come from PLINK PCA and Sex
    # case-control status (out_type = "D"),  "C" for the continuous outcome and "D" for the dichotomous outcome
    # SKAT_Null_Model will apply small sample size adjustment automatically when n<2000
    null_model <- SKAT_Null_Model(skat_data$phenotype_vector ~ 1, out_type = "D", Adjustment=TRUE)
    
    # method="optimal.adj" represent SKAT-O, default= "davies"
    # weights.beta use default value c(1,25)
    print("Start to run SKAT")
    set.seed(123)
    skat_results <- SKAT(skat_data$genotype_matrix, null_model,kernel="linear.weighted", method=kernel2use, weights.beta=c(1,25))
  }
  
  return(skat_results)
}

# -----------------------------------------------------------------------------
#### Usage ####
# -----------------------------------------------------------------------------
# Example Usage : Mutation burden of exonic and UTR variants in different geneList across DS-ECD and DS-nonECD samples
# Gene List to analysis from 04. Gene_in_genelists.R
# Prepare data
DS_ECD <- list(
  rbind(D25029_variantgroup$exonic, D25029_variantgroup$UTR),
  rbind(D25046_variantgroup$exonic, D25046_variantgroup$UTR)
)

DS_nonECD <- list(
  rbind(D25165_variantgroup$exonic, D25165_variantgroup$UTR),
  rbind(D25168_variantgroup$exonic, D25168_variantgroup$UTR)
)

SKATtest <- run_SKAT(gene_lists$top_candidate_related_genes, DS_ECD, DS_nonECD, "SKAT-O", by_gene = FALSE)
SKATtest$p.value
## FDR adjustment should be made if SKAT on each gene in genelist
adjusted_pvals <- p.adjust(pvals, method = "fdr")

## RUN SKAT BY GENE
# Run SKAT for each gene
SKATtest <- run_SKAT(gene_lists$top_candidate_related_genes, DS_ECD, DS_nonECD, "SKAT-O", by_gene = TRUE)

