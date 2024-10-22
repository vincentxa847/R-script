library(dplyr) 

#### Function to read tsv file that filter was applied using awk command ####
# Pre-Filter is Quality (Otherinfo9, col 267) ≥ 300 & Het (Het_Percent, col 359) > 15 & MAF (Max Allele Freq, 354) ≤0.01 | ==0 
ReadData <- function(filePath){
  
  # Generate the output name for the cleaned-up data
  dir <- dirname(filePath)
  file_name <- basename(filePath)
  file_prefix <- substr(file_name, 1, 6)  # Extract the first 6 characters, e.g., D25029
  output_name <- file.path(dir, paste0(file_prefix, ".MAF001.HET15.tsv"))  # Create output file path

  # Read the TSV file (allow missing values and non-standard quotes)
  # TSV file is generated using awk command, necessary to use [fill = TRUE, quote = ""]
  tmp <- read.table(filePath, header = TRUE, sep = "\t", fill = TRUE, quote = "")
  
  # Columns to move to the front
  front_columns <- c("Gene.refgene", "ExonicFunc.refgene", "CADD_phred", "REVEL_score", 
                     "Severity.Score", "Polyphen2_HDIV_score", "Polyphen2_HVAR_score",
                     "SIFT_score", "dbSNP", "ClinVar", "HGMD.Variant.Class", "Max.Allele.Freq", 
                     "Taiwan.Biobank", "Het_Percent", "Nucleotide", "AAChange")
  
  # Get all column names
  all_columns <- trimws(colnames(tmp))
  
  # Determine the other columns
  other_columns <- setdiff(all_columns, front_columns)
  
  # Reorder columns: front_columns first, followed by the remaining columns
  # Order by multiple columns
  tmp_combined_data <- tmp[, c(front_columns, other_columns)]
  combined_data <- tmp_combined_data[order(
    as.numeric(tmp_combined_data$CADD_phred), 
    as.numeric(tmp_combined_data$REVEL_score), 
    as.numeric(tmp_combined_data$Polyphen2_HDIV_score), 
    as.numeric(tmp_combined_data$SIFT_score),
    as.numeric(tmp_combined_data$SpliceAI_DS_AG),
    as.numeric(tmp_combined_data$SpliceAI_DS_AL),
    as.numeric(tmp_combined_data$SpliceAI_DS_DG),
    as.numeric(tmp_combined_data$SpliceAI_DS_DL),
    decreasing = TRUE
  ), ]
  
  # Variants in multiple overlapping genes were included in each gene separately 
  combined_data <- combined_data %>%
    tidyr::separate_rows(Gene.refgene, sep = ";")
  
  # Ambiguous annotation result in duplicated row, remove it
  combined_data <- combined_data[!duplicated(combined_data),]
  
  # Return the cleaned-up and reordered data (optional: write to output file)
  #write.table(combined_data, file = output_name, col.names=TRUE, row.names=FALSE, sep = "\t")
  combined_data
  
}

# DS-nonECD
D25163_MAF0.01 = ReadData("../WGS/HET15/D25163_22FHL2LT4.hg38_multianno.mviewer.tsv.MAF0.01.HET15.tsv")
D25165_MAF0.01 = ReadData("../WGS/HET15/D25165_22FHL2LT4.hg38_multianno.mviewer.tsv.MAF0.01.HET15.tsv")
D25168_MAF0.01 = ReadData("../WGS/HET15/D25168_22FHL2LT4.hg38_multianno.mviewer.tsv.MAF0.01.HET15.tsv")
# nonDS-ECD
D25007_MAF0.01 = ReadData("../WGS/HET15/D25007_22CYGGLT4.hg38_multianno.mviewer.tsv.MAF0.01.HET15.tsv")
# DS-ECD
D25029_MAF0.01 = ReadData("../WGS/HET15/D25029_22FHKYLT4.hg38_multianno.mviewer.tsv.MAF0.01.HET15.tsv") 
D25046_MAF0.01 = ReadData("../WGS/HET15/D25046_22FHKYLT4.hg38_multianno.mviewer.tsv.MAF0.01.HET15.tsv")
