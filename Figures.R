# Load necessary libraries
library(ggplot2)
library(pheatmap)
library(reshape2)
library(flextable)
library(officer)

# Set working directory
setwd("/Users/cassandradrazick/Desktop/6999Paper/Complete/")

# Read in data files for top 50 and top 2000 genes
top50 <- read.csv("Complete_50_gene_count_matrix.csv")
top2k <- read.csv("Complete_2k_gene_count_matrix.csv")

# Convert data to matrices
data_50 <- as.matrix(top50)
data_2k <- as.matrix(top2k)

# Set row names from the Gene column
rownames(data_50) <- top50$Gene
rownames(data_2k) <- top2k$Gene

# Remove the Gene column from the matrices
data_50 <- data_50[, -1]
data_2k <- data_2k[, -1]

# Convert matrices to data frames
data_2k <- as.data.frame(data_2k)

# Subset data_2k to include only genes with 'Total Occurrences' >= 10
subset_2k <- data_2k[data_2k$`Total Occurrences` >= 10, ]

# Ensure column names are consistent between datasets
colnames <- colnames(subset_2k)
colnames(data_50) <- colnames
colnames(subset_2k) <- colnames

# Reorder columns to place those with "C59" at the end
new_order <- c(
  setdiff(colnames, grep("C59", colnames, value = TRUE)),
  grep("C59", colnames, value = TRUE)
)
data_50 <- data_50[, new_order]
subset_2k <- subset_2k[, new_order]
subset_2k <- as.matrix(subset_2k)

## Generate heatmaps using ggplot2

# Melt data frames for ggplot compatibility
data_melt50 <- melt(data_50)
colnames(data_melt50) <- c("Gene", "Model", "Count")
data_melt50$Count <- as.numeric(as.character(data_melt50$Count))

data_melt2k <- melt(subset_2k)
colnames(data_melt2k) <- c("Gene", "Model", "Count")
data_melt2k$Count <- as.numeric(as.character(data_melt2k$Count))

# Define C59 column indices
unique_models <- unique(data_melt2k$Model)
c59_columns <- grep("C59", unique_models, value = TRUE)
c59_indices <- match(c59_columns, unique_models)

# Print C59 column indices for verification
print(c59_indices)

# Create heatmap for Top 50 genes
heatmap_50 <- ggplot(data_melt50, aes(x = Model, y = Gene, fill = Count)) +
  geom_tile() +
  scale_fill_viridis_c(option = "H") +
  theme_minimal() +
  labs(title = "Heatmap of Gene Patterns Occurring Greater Than Two Times, Top 50 Genes Across All Datasets", x = "Model", y = "Gene") +
  theme(axis.text.x = element_text(angle = 67, hjust = 1),
        axis.text.y = element_blank()) +
  geom_rect(data = data.frame(xmin = c59_indices - 0.5, xmax = c59_indices + 0.5,
                              ymin = -Inf, ymax = Inf),
            aes(x = NULL, y = NULL, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            color = "red", fill = NA, linewidth = 1)  

# Create heatmap for Top 2000 genes
heatmap_2k <- ggplot(data_melt2k, aes(x = Model, y = Gene, fill = Count)) +
  geom_tile() +
  scale_fill_viridis_c(option = "H") + 
  theme_minimal() +
  labs(title = "Heatmap of Gene Patterns Occurring Greater Than Ten Times, Top 2000 Genes Across All Datasets", x = "Model", y = "Gene") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank()) +
  geom_rect(data = data.frame(xmin = c59_indices - 0.5, xmax = c59_indices + 0.5,
                              ymin = -Inf, ymax = Inf),
            aes(x = NULL, y = NULL, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            color = "red", fill = NA, linewidth = 1)  

# Print the heatmaps
heatmap_50
heatmap_2k

### Create Flextables for Top 50 Genes across different datasets

# Set default font size for flextables
set_flextable_defaults(font.size = 6)

# Create and style flextables for different datasets
c59 <- flextable(C59_50_combined_top_genes_all_sources)
c59 <- set_caption(c59, caption = "Top 50 Genes Across Models for the C59 Dataset")
c59 <- set_table_properties(c59, width = 0.75, layout = "autofit")
c59 <- align(c59, align = "center", part = "all")

chir <- flextable(CHIR99021_combined_top_genes_all_sources)
chir <- set_caption(chir, caption = "Top 50 Genes Across Models for the CHIR99021 Dataset")
chir <- set_table_properties(chir, width = 0.75, layout = "autofit")
chir <- align(chir, align = "center", part = "all")

cebpa <- flextable(CEBPA_50_combined_top_genes_all_sources)
cebpa <- set_caption(cebpa, caption = "Top 50 Genes Across Models for the CEBPA Dataset")
cebpa <- set_table_properties(cebpa, width = 0.75, layout = "autofit")
cebpa <- align(cebpa, align = "center", part = "all")

fox <- flextable(FOXA1_50_combined_top_genes_all_sources)
fox <- set_caption(fox, caption = "Top 50 Genes Across Models for the FOXA1 Dataset")
fox <- set_table_properties(fox, width = 0.75, layout = "autofit")
fox <- align(fox, align = "center", part = "all")

axin <- flextable(nkdaxin_top50_genes_combined)
axin <- set_caption(axin, caption = "Top 50 Genes Across Models for the nkd1;axin2 Dataset")
axin <- set_table_properties(axin, width = 0.75, layout = "autofit")
axin <- align(axin, align = "center", part = "all")

WNT <- flextable(WNT3A_50_combined_top_genes_all_sources)
WNT <- set_caption(WNT, caption = "Top 50 Genes Across Models for the WNT3A Dataset")
WNT <- set_table_properties(WNT, width = 0.75, layout = "autofit")
WNT <- align(WNT, align = "center", part = "all")

# Create a Word document and add all flextables to it
doc1 <- read_docx()
doc1 <- body_add_flextable(doc1, c59)
doc1 <- body_add_flextable(doc1, chir)
doc1 <- body_add_flextable(doc1, cebpa)
doc1 <- body_add_flextable(doc1, axin)
doc1 <- body_add_flextable(doc1, fox)
doc1 <- body_add_flextable(doc1, WNT)
print(doc1, target = "Top_50_Genes1.docx")

### Create and style a flextable for subset of Top 50 Genes
data_50 <- as.data.frame(data_50)
data_50$`Total Occurrences` <- trimws(data_50$`Total Occurrences`)
data_50$`Total Occurrences` <- as.numeric(data_50$`Total Occurrences`)

# Check for any NA values or unexpected types in 'Total Occurrences'
print(str(data_50$`Total Occurrences`))
print(sum(is.na(data_50$`Total Occurrences`)))

# Subset data_50 with a threshold of 5
subset_50 <- data_50[data_50$`Total Occurrences` >= 5, ]

# Subset data_2k with a threshold of 16
subset_data2k <- as.data.frame(data_2k)
subset_data2k <- subset_data2k[subset_data2k$`Total Occurrences` >= 16, ]

# Create and style a flextable for occurrences of genes across models for Top 50 Genes
all <- flextable(subset_50)
all <- set_caption(all, caption = "Occurrences of Genes Across Models for Top 50 Genes")
all <- set_table_properties(all, width = 0.5, layout = "autofit")
all <- align(all, align = "center", part = "all")

### Create and style flextables for T-test results
setwd("/Users/cassandradrazick/Desktop/6999Paper/T-tests/")
C59_T <- read.delim("C59_50_ttests_results.csv", header = TRUE, stringsAsFactors = FALSE)
CHIR_T <- read.delim("CHIR99021_50_ttests_results.csv", header = TRUE, stringsAsFactors = FALSE)
CEBPA_T <- read.delim("CEBPA_50_ttests_results.csv", header = TRUE, stringsAsFactors = FALSE)
FOXA1_T <- read.delim("FOXA1_50_ttests_results.csv", header = TRUE, stringsAsFactors = FALSE)
NKD_T <- read.delim("nkd1;axin2_50_ttests_results.csv", header = TRUE, stringsAsFactors = FALSE)
WNT_T <- read.delim("WNT3A_50_ttests_results.csv", header = TRUE, stringsAsFactors = FALSE)

# Create flextables for T-test results and save them
C59_T <- flextable(C59_T)
C59_T <- set_caption(C59_T, caption = "T-test Results for Top 50 Genes, C59 Dataset")
C59_T <- set_table_properties(C59_T, width = 0.75, layout = "autofit")
C59_T <- align(C59_T, align = "center", part = "all")

CHIR_T <- flextable(CHIR_T)
CHIR_T <- set_caption(CHIR_T, caption = "T-test Results for Top 50 Genes, CHIR99021 Dataset")
CHIR_T <- set_table_properties(CHIR_T, width = 0.75, layout = "autofit")
CHIR_T <- align(CHIR_T, align = "center", part = "all")

CEBPA_T <- flextable(CEBPA_T)
CEBPA_T <- set_caption(CEBPA_T, caption = "T-test Results for Top 50 Genes, CEBPA Dataset")
CEBPA_T <- set_table_properties(CEBPA_T, width = 0.75, layout = "autofit")
CEBPA_T <- align(CEBPA_T, align = "center", part = "all")

FOXA1_T <- flextable(FOXA1_T)
FOXA1_T <- set_caption(FOXA1_T, caption = "T-test Results for Top 50 Genes, FOXA1 Dataset")
FOXA1_T <- set_table_properties(FOXA1_T, width = 0.75, layout = "autofit")
FOXA1_T <- align(FOXA1_T, align = "center", part = "all")

NKD_T <- flextable(NKD_T)
NKD_T <- set_caption(NKD_T, caption = "T-test Results for Top 50 Genes, nkd1;axin2 Dataset")
NKD_T <- set_table_properties(NKD_T, width = 0.75, layout = "autofit")
NKD_T <- align(NKD_T, align = "center", part = "all")

WNT_T <- flextable(WNT_T)
WNT_T <- set_caption(WNT_T, caption = "T-test Results for Top 50 Genes, WNT3A Dataset")
WNT_T <- set_table_properties(WNT_T, width = 0.75, layout = "autofit")
WNT_T <- align(WNT_T, align = "center", part = "all")

# Create a Word document and add all flextables to it
doc2 <- read_docx()
doc2 <- body_add_flextable(doc2, C59_T)
doc2 <- body_add_flextable(doc2, CHIR_T)
doc2 <- body_add_flextable(doc2, CEBPA_T)
doc2 <- body_add_flextable(doc2, FOXA1_T)
doc2 <- body_add_flextable(doc2, NKD_T)
doc2 <- body_add_flextable(doc2, WNT_T)
print(doc2, target = "Top_50_Genes_Ttests.docx")
