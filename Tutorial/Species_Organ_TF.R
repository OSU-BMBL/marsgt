# devtools::install_github("benoukraflab/TFregulomeR")
library(TFregulomeR)
library(DT)

## Function: Retrieve and save TF data
get_TF_data <- function(species, organ, output_file = "TFdata.csv") {
  # Load the TF database
  all_record <- dataBrowser()
  
  # Filter the database for the specified species and organ
  index <- which(all_record$species == species & all_record$organ == organ)
  
  # If no matching data is found, return an error message and stop
  if (length(index) == 0) {
    stop(paste("No data found for species:", species, "and organ:", organ))
  }
  
  # Extract relevant TF information
  TF_info <- all_record[index, c("ID", "species", "organ", "cell_tissue_name", "TF")]
  
  # Initialize a list to store TF Peaks data
  TFPeaks <- list()
  for (i in TF_info$ID) {
    TFPeaks[[i]] <- loadPeaks(id = i, includeMotifOnly = FALSE) # Load Peak data for each TF
  }
  
  # Combine all TF Peaks into a single data frame
  TF_dataframe <- do.call(rbind, TFPeaks)
  
  # Save the combined data frame as a CSV file
  write.table(TF_dataframe, file = output_file, row.names = FALSE, col.names = TRUE, sep = ",")
  
  # Return the filtered TF information and the combined Peaks data frame
  list(
    TF_info = TF_info,               # Filtered TF information
    TF_dataframe = TF_dataframe      # Combined TF Peaks data frame
  )
}

## Example: Call the function to retrieve data
result <- get_TF_data(species = "mouse", organ = "eye", output_file = "mouse_eye_TFdata.csv")

## Visualize the filtered TF information
DT::datatable(result$TF_info)
