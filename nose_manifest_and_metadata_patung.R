library(readr)
library(dplyr)

# Script to filter samples downloaded from HMP project and generate a manifest file to import fastaq sequences from HMP to QIIME2.

# Read metadata, includes sample_ids, samplee_body_location and visit_number to choose samples. Generated in "https://portal.hmpdacc.org/".
hmp_manifest_metadata_2df663e6ad <- read_delim("/srv/home/anavarro/hmp_nose/hmp_manifest_metadata_2df663e6ad.tsv", 
                                               delim = "\t", escape_double = FALSE, 
                                               trim_ws = TRUE)

# Read manifest file. This file is used to download fastq files from HMP servers using the portal_client package. Generated in "https://portal.hmpdacc.org/".
hmp_manifest_3d5b44c0d5 <- read_delim("/srv/home/anavarro/hmp_nose/hmp_manifest_3d5b44c0d5.tsv", 
                                      delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE)


# filter samples of "nasal cavity" and first subject's visit only in the metadata file.
# From this filtered df, we get the samples list and can be used to generate metadata df for these samples (for qiime2).
hmp_nose_metadata <- hmp_manifest_metadata_2df663e6ad %>% filter(sample_body_site == "nasal cavity")

# Get a vector of the samples id
sample_ids <- hmp_nose_metadata$sample_id

# Filter manifest file (where the names of the files are) to contain only the samples we chose previously.
filtered_manifest <- filter(hmp_manifest_3d5b44c0d5, sample_id %in% sample_ids)

# This function extracts the file name from the url in the manifest file. This file is a gz compressed file.
# GZ files are downloaded using portal_clinet and should be extracted in a separate directory each preserving the original name of the gz file.
get_file_name <- function(file_path_string){
  str_pieces <- strsplit(file_path_string, "/")[[1]]
  if (length(str_pieces) > 1) {
    file_name <- toString(str_pieces[length(str_pieces)])
    return(substr(file_name, 1, nchar(file_name)-3)) #remove the ".gz" part of the string
  }
  else{
    return(str_pieces[length(str_pieces)])
  }
}

# Get the names of the gz files that we are going to include in the qiime2 manifest file. Chosen using filters above.
files_names_to_include <- lapply(filtered_manifest$urls, get_file_name)

# Get the name of each directory result of the extraction of the gz files.
downloaded_files <- list.files("/srv/home/anavarro/hmp_nose/hmp_nose2", full.names = FALSE, recursive = FALSE)

# Create empty dataframe that is going to become the manifest file for qiime2
manifest_df <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(manifest_df) <- c("sample-id", "absolute-filepath")

# Iterate over all the sample ids that will be imported to qiime2
for(sample_number in seq(from = 1, to = length(sample_ids), by = 1)) {
  file_name <- files_names_to_include[sample_number]
  # If the sample id  is not already in the manifest file. This is because some samples have two or more fastaq files.
  if (!sample_ids[sample_number] %in% manifest_df$`sample-id` && !grepl(pattern = "Private", x = file_name, ignore.case = TRUE)) {
    if (file_name %in% downloaded_files) {
      absolute_file_path <- paste("/srv/home/anavarro/hmp_nose/hmp_nose2", file_name, sep = "/") # get string with directory and fastaq file name for each sample.
      manifest_df[nrow(manifest_df) + 1,] = c(sample_ids[sample_number], absolute_file_path)
    }
  }
}


write.table(manifest_df, file = "/srv/home/anavarro/hmp_nose/manifest_nasal_cavity_all.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
