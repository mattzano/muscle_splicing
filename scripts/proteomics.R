###proteomics search
library(data.table)
library(tidyverse)

# Set the directory containing the files
directory <- here::here('data',"68k_entries_sample-FDR/de novo")

# List all files in the directory containing "de_novo" in their names
file_list <- list.files(path = directory, pattern = "AllCandidates", full.names = TRUE)

# Function to read each CSV file
read_csv_file <- function(file) {
  read.csv(file, stringsAsFactors = FALSE)
}

# Read all the files into a list of dataframes
dataframes <- lapply(file_list, read_csv_file)

# Combine all dataframes into one
combined_dataframe <- do.call(rbind, dataframes)

# View the combined dataframe
#print(combined_dataframe)


peptides <- combined_dataframe$Peptide

# String to search for
search_string <- ""
#EKKAAVRAPRRGPLGGRKKK EPTIWFGKGHSGMLASEGREAVLTRLHESERVRKQERERDTEERRE KAPSASDSDSKADSDGAKPE"

# Using grep to get indices of peptides that contain the search string
matching_indices <- grep(search_string, peptides)

# Using grepl to get a logical vector indicating presence of the search string
matching_logical <- grepl(search_string, peptides)
matching_peptides <- peptides[matching_logical]




prot_table <- read.csv(here::here("data", "68k_entries_sample-FDR/dia_db.proteins.csv"))

prot_table %>% 
  filter(Gene == "TARDBP") %>% 
  pivot_longer(cols = c(7:27)) %>% 
  ggplot(aes(x = reorder(name, value), y = value)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))



data <- read.csv("~/Downloads/dia_db.peptides_Cryptic.csv")

# Add '>' to accession values
data$Accession <- paste0(">", data$Accession)

# Reshape the data into a single column
long_data <- pivot_longer(data, cols = c(2:1))

fastino <- long_data$value

write.table(fastino, "~/Desktop/ce_pep.fa", quote = F, row.names = F)





# Define the two lists of strings
list1 <- read.csv("~/Downloads/manual_peptides.xlsx - Sheet 1 - manual_peptides_57_.csv") %>% 
  pull(flanked)
list2 <- read.csv("~/Downloads/dia_db.peptides_Cryptic.csv") %>% 
  pull(Peptide)
#list1 <- c("apple", "banana", "cherry")
#list2 <- c("app", "nan", "err", "fruit")

# Create an empty matrix to store match results
match_matrix <- matrix(FALSE, nrow = length(list1), ncol = length(list2))

# Check for partial matches
for (i in seq_along(list1)) {
  for (j in seq_along(list2)) {
    match_matrix[i, j] <- grepl(list2[j], list1[i])
  }
}

# Convert the matrix to a data frame for better readability
match_df <- as.data.frame(match_matrix, row.names = list1, col.names = list2)

