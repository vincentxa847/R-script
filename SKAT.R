library(SKAT)

#### Prepare the Data for SKAT (Genotype Matrix, Phenotype Vector) #### 
## Function to extract gene in gene list and corresponding genotype from variant table 
extract_genotype_filtered <- function(sample_data, gene_list) {
  # Filter the sample data to only include rows where Gene_refgene matches the gene list
  filtered_data <- sample_data[sample_data$Gene_refgene %in% gene_list, ]
  
  # Create a variant identifier by combining "Chr" and "Start"
  filtered_data$Variant <- paste(filtered_data$Chr, filtered_data$Start, sep = "_")
  
  # Extract genotype and convert 'het' -> 1 and 'hom' -> 2
  # Assuming homozygous for the minor allele (two copies of the variant allele) is represented as 2
  genotypes <- ifelse(filtered_data$Otherinfo1 == "het", 1, 2)
  
  # Return a data frame with the new variant identifier and genotypes
  return(data.frame(Variant = filtered_data$Variant, Genotype = genotypes))
}

## Function to create genotype matrix and phenotype vector for SKAT
generate_SKAT_data <- function(case_list, control_list, gene_list) {
  
  # Create a vector of genes from the gene list table
  gene_list <- gene_list$Gene_refgene
  
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

#### SKAT Analysis of the mutation burden of rare non-synonymous variants (missense) in cilium components across DS-ECD and non-DS-ECD samples ####
## Gene List to analysis (cilium_components in this case)
cilium_components <- read.csv("../GENE_LIST/cilium_components.txt",
                              header = FALSE, sep = "\t", 
                              stringsAsFactors = FALSE, col.names = "Gene_refgene")

## Combine rare non-synonymous variants of samples of each group into a list 
DS_ECD <- list(D25029_MAF0.01_missense, D25046_MAF0.01_missense)
nonDS_ECD <- list(D25007_MAF0.01_missense)

##  Define covariate (for now sample size is too small (n=3) to include covariate)
# include sex as covariate (1 = Male, 0 = Female)
sex_vector <- c(1, 0, 0)  # D25029:M, D25046:F, D25007:F
# Create a data frame for covariates
covariate_df <- data.frame(sex = sex_vector)

## Generate the genotype matrix and phenotype vector for SKAT analysis
skat_data <- generate_SKAT_data(DS_ECD, nonDS_ECD, cilium_components)

## Run SKAT-O
# Creating the Null Model, which represents the baseline against which the alternative hypothesis
# Not using covariate here (~ 1), which come from PLINK PCA and Sex
# case-control status (out_type = "D"),  "C" for the continuous outcome and "D" for the dichotomous outcome
null_model <- SKAT_Null_Model(skat_data$phenotype_vector ~ 1, out_type = "D")
# method="optimal.adj" represent SKAT-O, default= "davies"
# weights.beta use default value c(1,25)
skat_result <- SKAT(skat_data$genotype_matrix, null_model,kernel="linear.weighted", method="optimal.adj", weights.beta=c(1,25))
