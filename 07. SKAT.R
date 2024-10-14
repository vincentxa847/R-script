library(SKAT)

# -----------------------------------------------------------------------------
#### Prepare the Data for SKAT (Genotype Matrix, Phenotype Vector) #### 
# -----------------------------------------------------------------------------
## Function to extract gene in gene list and corresponding genotype from variant table 
extract_genotype_filtered <- function(sample_data, gene_list) {
  # Filter the sample data to only include rows where Gene_refgene matches the gene list
  filtered_data <- sample_data[sample_data$Gene_refgene %in% gene_list, ]
  
  # Create a variant identifier by combining "Chr" and "Start"
  filtered_data$Variant <- paste(filtered_data$Chr, filtered_data$Start, sep = "_")
  
  # Extract genotype and convert 'het' -> 1 and 'hom' -> 2
  # Assuming homozygous for the minor allele (two copies of the variant allele) is represented as 2
  # A suggested approach is to set "hem" (hemizygous) to "homo" (2) or exclude (9). For now, it's set to 2.
  genotypes <- ifelse(filtered_data$Otherinfo1 == "het", 1, 2)
  
  # Return a data frame with the new variant identifier and genotypes
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
    sample_genotype <- extract_genotype_filtered(case_list[[i]], gene_list)
    genotype_list[[paste0("Case_", i)]] <- sample_genotype
    all_variants <- union(all_variants, sample_genotype$Variant)  
    sample_names <- c(sample_names, paste0("Case_", i))
  }
  
  # Process control data
  for (i in seq_along(control_list)) {
    sample_genotype <- extract_genotype_filtered(control_list[[i]], gene_list)
    genotype_list[[paste0("Control_", i)]] <- sample_genotype
    all_variants <- union(all_variants, sample_genotype$Variant)  
    sample_names <- c(sample_names, paste0("Control_", i))
  }
  
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
  
  # Return the genotype matrix and phenotype vector
  return(list(genotype_matrix = genotype_matrix, phenotype_vector = phenotype_vector))
}

## Function to run SKAT
run_SKAT <- function(geneList, group1, group2, SKAT_or_SKATO){
  
  genelist2analysis <- geneList
  
  if (SKAT_or_SKATO == "SKAT") {
    method2use <- "davies"  # SKAT method
  } else if (SKAT_or_SKATO == "SKAT-O") {
    method2use <- "optimal.adj"  # SKAT-O method
  }
  
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
  set.seed(123)
  skat_result <- SKAT(skat_data$genotype_matrix, null_model, kernel="linear.weighted", method=kernel2use, weights.beta=c(1,25))
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

SKATtest <- run_SKAT(gene_lists$cilium, DS_ECD, DS_nonECD, "SKAT-O")
