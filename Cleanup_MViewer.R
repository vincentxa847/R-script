## Clean up data from MViewer  
MViewer_cleanup <-function(input_filePath){
  
  output_name <- sub("0.01", "0.01-cleanup", input_filePath)
  
  # Read the file line by line 
  lines <- readLines(input_filePath)
  # Initialize an empty list to store processed rows
  processed_data <- list()
  
  # Loop over each line
  for (i in seq_along(lines)) {
    # Split the line into columns based on tab delimiter
    split_line <- strsplit(lines[i], "\t")[[1]]
    
    # Keep only the first 372 columns, or fewer if the line has less
    split_line <- split_line[1:min(372, length(split_line))]
    
    # Add the processed line to the list
    processed_data[[i]] <- split_line
  }
  
  # Combine the list into a data frame
  combined_data <- as.data.frame(do.call(rbind, processed_data), stringsAsFactors = FALSE)
  
  # Set the first row as the column names
  colnames(combined_data) <- combined_data[1, ]
  
  # Remove the first row, which is now the column names
  combined_data <- combined_data[-1, ]
  
  
  # Ensure column names are trimmed and then check for the presence of the required columns
  if (all(c("Gene_refgene", "ExonicFunc_refgene", "CADD_phred", "REVEL_score") %in% trimws(colnames(combined_data)))) {
    # Define the columns to move to the front
    front_columns <- c("Gene_refgene", "ExonicFunc_refgene", "CADD_phred", "REVEL_score", "Severity_Score", "dbSNP", "ClinVar", "Max_Allele_Freq", "Taiwan_Biobank","Nucleotide" ,"AAChange")
    
    # Get all column names
    all_columns <- trimws(colnames(combined_data))
    
    # Define the columns to keep (excluding "Item" and columns already in front_columns)
    other_columns <- setdiff(all_columns, c(front_columns, "Item"))
    
    # Reorder columns
    combined_data <- combined_data[, c(front_columns, other_columns)]
  }
  combined_data <- combined_data[order(as.numeric(combined_data$CADD_phred), decreasing = TRUE), ]
  
  
  # Overwrite the original data to cleanup data
  write.table(combined_data, file = output_name, col.names=TRUE, row.names=FALSE, sep = "\t")
  
  combined_data
}

# OUTPUT THE ORIGINAL TABLE FROM MVIEWER (ONLY MAF<0.01) 
N1675_MAF0.01 = MViewer_cleanup("../N1675/N1675-MAF0.01.tsv")
AI3008_MAF0.01 = MViewer_cleanup("../AI3008/AI3008-MAF0.01.tsv")
# OUTPUT THE ORIGINAL TABLE FROM MVIEWER (ONLY MAF<0.01) 
N1675_MAF0.01 = MViewer_cleanup("D:/DS Project_Vincent/Data/N1675/N1675-MAF0.01.tsv")
AI3008_MAF0.01 = MViewer_cleanup("D:/DS Project_Vincent/Data/AI3008/AI3008-MAF0.01.tsv")
