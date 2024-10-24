library(dplyr)


# -------------------------------------------------------------------
# Usage of expression data
# Reference :  
# -------------------------------------------------------------------

# Step 1: List all files in the directory
files <- list.files(path = "../Expression_Data/GSE106118_RAW/", pattern = "AV.*\\.txt\\.gz$", full.names = TRUE)

# Step 2: Filter files that have the specific time points in the filename
time_points <- c("5W", "6W", "7W", "9W", "10W", "13W", "15W", "17W", "20W", "22W", "23W", "24W", "25W")
pattern <- paste0(time_points, collapse = "|") # Create regex pattern with time points

# four expression file for AV
#[1] "../Expression_Data/GSE106118_RAW/GSM2829989_HE22W_AV.TPM.txt.gz"    "../Expression_Data/GSE106118_RAW/GSM2830005_HE23W_2_AV.TPM.txt.gz" 
#[3] "../Expression_Data/GSE106118_RAW/GSM2830015_HE25W_AV_EP.TPM.txt.gz" "../Expression_Data/GSE106118_RAW/GSM3208696_HE17W_3_AV_TPM.txt.gz" 
files_filtered <- files[grepl(pattern, files)]

# Step 3: Function to read and process each file
process_file <- function(file) {
  # Read in the gzipped file (assuming tab-delimited)
  data <- read.table(gzfile(file), header = TRUE, row.names = 1)  # Adjust delimiter if necessary
  
  # Step 4: Sum the rows of the data
  row_sums <- rowSums(data) # or you can use rowsum if grouping is needed
  
  return(row_sums)
}

# Step 4: Apply the process to each file and sum the results
# This reduces (combines) the list of vectors by summing them element-wise. In this case, it sums the row sums of all files into a single vector
total_sum <- Reduce(`+`, lapply(files_filtered, process_file))

# Access expression value
total_sum["GRIN3B"]
