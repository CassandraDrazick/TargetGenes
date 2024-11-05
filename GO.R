# Load necessary libraries
library(dbplyr)
library(dplyr)
library(biomaRt)
library(gprofiler2)
library(wds)  

##CHIR99021

# Load data files
norm_counts_no_zeros <- read.table("norm_counts_no_zeros.tsv", header = TRUE, sep = "\t", row.names = 1)
gene_names <- read.table("gene_names.tsv", header = TRUE, sep = "\t")

# Step 1: Remove decimals from gene IDs in the norm_counts_no_zeros dataframe
gene_ids <- gsub("\\..*$", "", row.names(norm_counts_no_zeros))

# Step 2: Prepare gene_names_df to match the updated gene_ids
gene_names_df <- data.frame(
  gene_id = gene_names$input,
  gene_name = gene_names$name
)

# Ensure the gene_id in gene_names_df is a character vector
gene_names_df$gene_id <- as.character(gene_names_df$gene_id)

# Step 3: Identify gene IDs in gene_ids that are not in gene_names_df
missing_genes <- setdiff(gene_ids, gene_names_df$gene_id)

# Step 4: Remove missing gene IDs from gene_ids and norm_counts_no_zeros
filtered_gene_ids <- gene_ids[!gene_ids %in% missing_genes]
filtered_norm_counts_no_zeros <- norm_counts_no_zeros[rownames(norm_counts_no_zeros) %in% filtered_gene_ids, ]

# Step 5: Merge filtered_norm_counts_no_zeros with gene_names_df
CCantu_geneIDs <- merge(
  filtered_norm_counts_no_zeros,
  gene_names_df,
  by.x = "row.names",
  by.y = "gene_id",
  all.x = TRUE
)

# Step 6: Assign rownames and clean up CCantu_geneIDs
rownames(CCantu_geneIDs) <- CCantu_geneIDs$`Row.names`
CCantu_geneIDs$`Row.names` <- NULL

# Step 7: Identify duplicates and combine them
CCantu_geneIDs_combined <- CCantu_geneIDs %>%
  group_by(gene_name) %>%
  summarise(across(everything(), mean))

# Convert to a regular dataframe
CCantu_geneIDs_combined <- as.data.frame(CCantu_geneIDs_combined)

# Set gene_name as the row names
rownames(CCantu_geneIDs_combined) <- CCantu_geneIDs_combined$gene_name

# Remove the gene_name column
CCantu_geneIDs_combined$gene_name <- NULL

# Step 8: View the updated dataframe
View(CCantu_geneIDs_combined)  # View dataframe in RStudio or print to console

# Step 9: Save the final dataframe to a file
write.table(CCantu_geneIDs_combined, file = "finalCCantu.tsv", sep = "\t", row.names = TRUE, col.names = NA)

##Nkd1Axin2

# Connect to the Ensembl database and specify the dataset
ensembl <- useMart("ensembl", dataset = "drerio_gene_ensembl")
datasets <- listDatasets(ensembl)
head(datasets)

# Retrieve gene mappings for zebrafish ensembl_gene_id and external_gene_name
gene_mappings <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), values = zebrafish_genes, mart = ensembl)

# Load normalized data file
normalized_cleaned_IAN <- read.table("normalized_cleaned_IAN.tsv", header = TRUE, sep = "\t", row.names = 1)

# Add rownames as a new column
normalized_cleaned_IAN$ensembl_gene_id <- rownames(normalized_cleaned_IAN)

# Merge the dataframes
merged_df <- merge(normalized_cleaned_IAN, gene_mappings, by = "ensembl_gene_id")

# Reorder columns if needed
merged_df <- merged_df[, c("human_gene_name", setdiff(colnames(merged_df), c("human_gene_name", "ensembl_gene_id")))]

# Remove the ensembl_gene_id column if no longer needed
merged_df <- merged_df[, !colnames(merged_df) %in% "ensembl_gene_id"]

# Merge gene_mappings with human_gene_names
mapping_table <- merge(
  gene_mappings,
  human_gene_names,
  by.x = "hsapiens_homolog_ensembl_gene",
  by.y = "ensembl_gene_id"
)

# Create a lookup table for converting zebrafish ENS IDs to human gene names
lookup_table <- mapping_table[, c("ensembl_gene_id", "external_gene_name.y")]
colnames(lookup_table) <- c("zebrafish_ensembl_gene_id", "human_gene_name")

# Match zebrafish ENS IDs to human gene names
new_rownames <- lookup_table$human_gene_name[match(rownames(normalized_cleaned_IAN), lookup_table$zebrafish_ensembl_gene_id)]

# Add new_rownames as a column in the normalized_cleaned_IAN dataframe
df <- as.data.frame(normalized_cleaned_IAN)
df$human_gene_name <- new_rownames

# Check for and handle missing values
missing_genes <- which(is.na(df$human_gene_name))
if (length(missing_genes) > 0) {
  warning("Some zebrafish ENS IDs do not have corresponding human gene names: ", paste(rownames(normalized_cleaned_IAN)[missing_genes], collapse = ", "))
}

# Remove rows with NA human_gene_name
df <- df[!is.na(df$human_gene_name), ]

# Aggregate the data by calculating the mean for duplicate gene names
df_aggregated <- df %>%
  group_by(human_gene_name) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
df_aggregated <- as.data.frame(df_aggregated)
# Convert the aggregated data frame back to a matrix
df_aggregated <- df_aggregated %>%
  column_to_rownames(var = "human_gene_name") 

# Reset row names if necessary
rownames(df_aggregated) <- df_aggregated$human_gene_name
df_aggregated <- df_aggregated[, -which(colnames(df_aggregated) == "human_gene_name")]

# Write the human_genes_IAN data frame to a TSV file
write.table(human_genes_IAN, file = "human_genes_IAN.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

##WNT3A and C59

# Define a function to convert and aggregate gene codes
convert_and_aggregate <- function(df) {
  gene_codes <- rownames(df)
  
  # Use gconvert function from gprofiler2 to map gene codes to gene names
  gene_annotations <- gconvert(query = gene_codes, organism = "hsapiens", target = "ENSG", filter_na = TRUE)
  
  if (nrow(gene_annotations) == 0) {
    stop("No gene annotations were retrieved. Please check the gene codes and try again.")
  }
  
  print("Gene annotations retrieved:")
  print(head(gene_annotations))
  
  # Merge the gene names with the original dataframe
  df$Gene <- rownames(df)
  annotated_df <- merge(df, gene_annotations, by.x = "Gene", by.y = "input", all.x = TRUE)
  
  if (nrow(annotated_df) == 0) {
    stop("The merge operation resulted in an empty dataframe. Please check the merge step.")
  }
  
  print("Merged dataframe:")
  print(head(annotated_df))
  
  # Ensure numeric columns are correctly identified for summarization
  numeric_cols <- sapply(annotated_df, is.numeric)
  
  # Aggregate duplicates by taking the mean of the expression values
  aggregated_df <- annotated_df %>%
    group_by(name) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = 'drop')
  
  if (nrow(aggregated_df) == 0) {
    stop("The aggregation step resulted in an empty dataframe. Please check the aggregation step.")
  }
  
  print("Aggregated dataframe:")
  print(head(aggregated_df))
  
  # Set row names to gene names
  aggregated_df <- aggregated_df[!is.na(aggregated_df$name), ]
  aggregated_df <- as.data.frame(aggregated_df) # Convert to dataframe
  rownames(aggregated_df) <- aggregated_df$name
  aggregated_df <- aggregated_df[, -which(names(aggregated_df) %in% c("name", "Gene"))]
  
  if (nrow(aggregated_df) == 0 || ncol(aggregated_df) == 0) {
    stop("The final dataframe is empty. Please check the row name setting and column removal steps.")
  }
  
  return(aggregated_df)
}

# Load data files
combined_ES_WNT <- read.table("combined_ES_WNT.tsv", header = TRUE, sep = "\t", row.names = 1)
combined_ES_C59 <- read.table("combined_ES_C59.tsv", header = TRUE, sep = "\t", row.names = 1)

# Convert gene codes and aggregate duplicates for combined_ES_WNT
combined_ES_WNTgene <- convert_and_aggregate(combined_ES_WNT)

# Convert gene codes and aggregate duplicates for combined_ES_C59
combined_ES_C59gene <- convert_and_aggregate(combined_ES_C59)

# Save cleaned dataframes
write.table(combined_ES_WNTgene, file = "combined_ES_WNT_clean.tsv", sep = "\t", row.names = TRUE, quote = FALSE)
write.table(combined_ES_C59gene, file = "combined_ES_C59_clean.tsv", sep = "\t", row.names = TRUE, quote = FALSE)
