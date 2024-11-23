# Load required libraries
# install.packages("qs")
library(qs)
library(GenomicRanges)

# Define a function to check overlaps between peaks and TF binding sites
findTFOverlaps <- function(jaspar_path, file_name, peaks) {
  # Construct the full file path
  file_path <- paste0(jaspar_path, file_name)
  
  # Check if the file exists
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  # Load the .qsave file containing TF binding site data
  cat("Loading qsave file...\n")
  tfbs_df <- qs::qread(file_path)
  
  # Assign column names to the data frame
  colnames(tfbs_df) <- c("chr", "start", "end", "TF", "score", "p_value", "strand")
  
  # Convert the peaks data frame to a GRanges object
  peak_gr <- makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
  
  # Convert the TF binding sites data frame to a GRanges object
  tf_binding_gr <- makeGRangesFromDataFrame(tfbs_df, keep.extra.columns = TRUE, na.rm = TRUE)
  
  # Use findOverlaps to determine overlapping regions
  overlaps <- findOverlaps(query = peak_gr, subject = tf_binding_gr)
  
  # Extract the overlapping TF binding sites
  overlap_results <- tf_binding_gr[subjectHits(overlaps)]
  
  # Return the overlapping results
  return(overlap_results)
}

# Example usage:
# Define the path to the JASPAR data
jaspar_path <- "jaspar_data/"
file_name <- "hg38_lisa_500.qsave" # Human data (hg38 genome)
# file_name <- "mm10.qsave"        # Mouse data (mm10 genome)

# Define the peaks of interest (user-customized)
peaks <- data.frame(
  chr = "chr10",
  start = 28000,
  end = 28600
)

# Call the function to find overlaps
overlap_results <- findTFOverlaps(jaspar_path, file_name, peaks)

# Print the results
print(overlap_results)
