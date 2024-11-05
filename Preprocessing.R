#####Preprocessing Steps for each dataset#####

#####C59 and WNT3A Muckerjee et al.(2022)#####

setwd("/Users/cassandradrazick/Desktop/BINF6999/Muckerjee")

####FASTQ

# Load packages
library(Rsubread)
library(edgeR)
library(tidyverse)  
library(biomaRt)
library(dbplyr)
library(readr)

##Preprocess
##Feature counts 
# Define the annotation file (GTF) THE FILE WAS SO MASSIVE I HAD TO DO IT SEPARATELY
annotation_file <- "/Users/cassandradrazick/Downloads/Muckerjee/Homo_sapiens.GRCh38.112.gtf.gz"
es_bam_file1 <-  "SRR15622185_Aligned.sortedByCoord.out.bam"
es_bam_file2 <- "SRR15622186_Aligned.sortedByCoord.out.bam"
es_bam_file3 <-"SRR15622187_Aligned.sortedByCoord.out.bam"
wnt1_bam_file1 <- "SRR15622188_Aligned.sortedByCoord.out.bam"
wnt1_bam_file2 <- "SRR15622189_Aligned.sortedByCoord.out.bam"
wnt1_bamfile3 <- "SRR15622190_Aligned.sortedByCoord.out.bam"
wnt2_bam_file1 <- "SRR15622194_Aligned.sortedByCoord.out.bam"
wnt2_bam_file2 <- "SRR15622195_Aligned.sortedByCoord.out.bam"
wnt2_bam_file3 <- "SRR15622196_Aligned.sortedByCoord.out.bam"
c592e_bam_file1 <- "SRR15622197_Aligned.sortedByCoord.out.bam"
c592e_bam_file2 <- "SRR15622198_Aligned.sortedByCoord.out.bam"
c592e_bam_file3 <- "SRR15622199_Aligned.sortedByCoord.out.bam"
c592l_bam_file1 <- "SRR15622200_Aligned.sortedByCoord.out.bam"
c592l_bam_file2 <- "SRR15622201_Aligned.sortedByCoord.out.bam"
c592l_bam_file3 <- "SRR15622202_Aligned.sortedByCoord.out.bam"
wnt3_bam_file1 <- "SRR15622203_Aligned.sortedByCoord.out.bam"
wnt3_bam_file2 <- "SRR15622204_Aligned.sortedByCoord.out.bam"
wnt3_bam_file3 <- "SRR15622205_Aligned.sortedByCoord.out.bam"
c593e_bam_file1 <- "SRR15622206_Aligned.sortedByCoord.out.bam"
c593e_bam_file2 <- "SRR15622207_Aligned.sortedByCoord.out.bam"
c593e_bam_file3 <- "SRR15622208_Aligned.sortedByCoord.out.bam"
c593l_bam_file1 <- "SRR15622209_Aligned.sortedByCoord.out.bam"
c593l_bam_file2 <- "SRR15622210_Aligned.sortedByCoord.out.bam"
c593l_bam_file3 <- "SRR15622211_Aligned.sortedByCoord.out.bam"

#Count features to get gene expression levels
fcMUCK_ES1 <- featureCounts(es_bam_file1, annot.ext = annotation_file, isGTFAnnotationFile = TRUE, 
                            GTF.featureType = "exon", GTF.attrType = "gene_id", 
                            useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)
fcMUCK_ES2 <- featureCounts(es_bam_file2, annot.ext = annotation_file, isGTFAnnotationFile = TRUE, 
                            GTF.featureType = "exon", GTF.attrType = "gene_id", 
                            useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)
fcMUCK_ES3 <- featureCounts(es_bam_file3, annot.ext = annotation_file, isGTFAnnotationFile = TRUE, 
                            GTF.featureType = "exon", GTF.attrType = "gene_id", 
                            useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)


# Extract the counts from each featureCounts output object
counts_ES1 <- fcMUCK_ES1$counts
counts_ES2 <- fcMUCK_ES2$counts
counts_ES3 <- fcMUCK_ES3$counts

# Combine the counts into a single dataframe
combined_counts <- cbind(counts_ES1, counts_ES2, counts_ES3)

# Rename the columns for clarity
colnames(combined_counts) <- c("ES1", "ES2", "ES3")

# Set the rownames to gene names or IDs if they are not already set
rownames(combined_counts) <- rownames(fcMUCK_ES1$counts)  # Uncomment and adjust if needed


##WNT feature counts
fcMUCK_WNT11 <- featureCounts(wnt1_bam_file1, annot.ext = annotation_file, isGTFAnnotationFile = TRUE, 
                              GTF.featureType = "exon", GTF.attrType = "gene_id", 
                              useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)
fcMUCK_WNT12 <- featureCounts(wnt1_bam_file2, annot.ext = annotation_file, isGTFAnnotationFile = TRUE, 
                              GTF.featureType = "exon", GTF.attrType = "gene_id", 
                              useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)
fcMUCK_WNT13 <- featureCounts(wnt1_bamfile3, annot.ext = annotation_file, isGTFAnnotationFile = TRUE, 
                              GTF.featureType = "exon", GTF.attrType = "gene_id", 
                              useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)

# Extract the counts from each featureCounts output object
counts_WNT11 <- fcMUCK_WNT11$counts
counts_WNT12 <- fcMUCK_WNT12$counts
counts_WNT13 <- fcMUCK_WNT13$counts

# Combine the counts into a single dataframe
combined_counts_WNT <- cbind(counts_WNT11, counts_WNT12, counts_WNT13)

# Rename the columns for clarity
colnames(combined_counts_WNT) <- c("WNT1_1", "WNT1_2", "WNT1_3")

# Optionally, set the rownames to gene names or IDs if they are not already set
rownames(combined_counts_WNT) <- rownames(fcMUCK_WNT11$counts)  

# View the combined counts dataframe
print(head(combined_counts_WNT))

##WNT2
fcMUCK_WNT21 <- featureCounts(wnt2_bam_file1, annot.ext = annotation_file, isGTFAnnotationFile = TRUE, 
                              GTF.featureType = "exon", GTF.attrType = "gene_id", 
                              useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)
fcMUCK_WNT22 <- featureCounts(wnt2_bam_file2, annot.ext = annotation_file, isGTFAnnotationFile = TRUE, 
                              GTF.featureType = "exon", GTF.attrType = "gene_id", 
                              useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)
fcMUCK_WNT23 <- featureCounts(wnt2_bam_file3, annot.ext = annotation_file, isGTFAnnotationFile = TRUE, 
                              GTF.featureType = "exon", GTF.attrType = "gene_id", 
                              useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)

# Extract the counts from each featureCounts output object
counts_WNT21 <- fcMUCK_WNT21$counts
counts_WNT22 <- fcMUCK_WNT22$counts
counts_WNT23 <- fcMUCK_WNT23$counts

# Combine the counts into a single dataframe
combined_counts_WNT2 <- cbind(counts_WNT21, counts_WNT22, counts_WNT23)

# Rename the columns for clarity
colnames(combined_counts_WNT2) <- c("WNT2_1", "WNT2_2", "WNT2_3")

# Set the rownames to gene names or IDs if they are not already set
rownames(combined_counts_WNT2) <- rownames(fcMUCK_WNT21$counts)  

##WNT3

fcMUCK_WNT31 <- featureCounts("SRR15622203_Aligned.sortedByCoord.out.bam", annot.ext = annotation_file, isGTFAnnotationFile = TRUE, 
                              GTF.featureType = "exon", GTF.attrType = "gene_id", 
                              useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)
fcMUCK_WNT32 <- featureCounts("SRR15622204_Aligned.sortedByCoord.out.bam", annot.ext = annotation_file, isGTFAnnotationFile = TRUE, 
                              GTF.featureType = "exon", GTF.attrType = "gene_id", 
                              useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)

fcMUCK_WNT33 <- featureCounts("SRR15622205_Aligned.sortedByCoord.out.bam", annot.ext = annotation_file, isGTFAnnotationFile = TRUE, 
                              GTF.featureType = "exon", GTF.attrType = "gene_id", 
                              useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)

# Extract the counts from each featureCounts output object
counts_WNT31 <- fcMUCK_WNT31$counts
counts_WNT32 <- fcMUCK_WNT32$counts
counts_WNT33 <- fcMUCK_WNT33$counts

# Combine the counts into a single dataframe
combined_counts_WNT3 <- cbind(counts_WNT31, counts_WNT32, counts_WNT33)

# Rename the columns for clarity
colnames(combined_counts_WNT3) <- c("WNT3_1", "WNT3_2", "WNT3_3")

# Set the rownames to gene names or IDs if they are not already set
rownames(combined_counts_WNT3) <- rownames(fcMUCK_WNT31$counts)  

##C591
fcMUCK_C592E1 <- featureCounts(c592e_bam_file1, annot.ext = annotation_file, isGTFAnnotationFile = TRUE,
                               GTF.featureType = "exon", GTF.attrType = "gene_id", 
                               useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)
fcMUCK_C592E2 <- featureCounts(c592e_bam_file2, annot.ext = annotation_file, isGTFAnnotationFile = TRUE,
                               GTF.featureType = "exon", GTF.attrType = "gene_id", 
                               useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)
fcMUCK_C592E3 <- featureCounts(c592e_bam_file3, annot.ext = annotation_file, isGTFAnnotationFile = TRUE,
                               GTF.featureType = "exon", GTF.attrType = "gene_id", 
                               useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)

# Extract the counts from each featureCounts output object
counts_C592E1 <- fcMUCK_C592E1$counts
counts_C592E2 <- fcMUCK_C592E2$counts
counts_C592E3 <- fcMUCK_C592E3$counts

# Combine the counts into a single dataframe
combined_counts_WNT3 <- cbind(counts_WNT31, counts_WNT32, counts_WNT33)

# Rename the columns for clarity
colnames(combined_counts_WNT3) <- c("WNT3_1", "WNT3_2", "WNT3_3")

# Set the rownames to gene names or IDs if they are not already set
rownames(combined_counts_WNT3) <- rownames(fcMUCK_WNT31$counts)  

##C592
fcMUCK_C592L1 <- featureCounts(c592l_bam_file1, annot.ext = annotation_file, isGTFAnnotationFile = TRUE,
                               GTF.featureType = "exon", GTF.attrType = "gene_id", 
                               useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)
fcMUCK_C592L2 <- featureCounts(c592l_bam_file2, annot.ext = annotation_file, isGTFAnnotationFile = TRUE,
                               GTF.featureType = "exon", GTF.attrType = "gene_id", 
                               useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)
fcMUCK_C592L3 <- featureCounts(c592l_bam_file3, annot.ext = annotation_file, isGTFAnnotationFile = TRUE,
                               GTF.featureType = "exon", GTF.attrType = "gene_id", 
                               useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)
# Extract the counts from each featureCounts output object
counts_C592L1 <- fcMUCK_C592L1$counts
counts_C592L2 <- fcMUCK_C592L2$counts
counts_C592L3 <- fcMUCK_C592L3$counts

# Combine the counts into a single dataframe
combined_counts_C592L <- cbind(counts_C592L1, counts_C592L2, counts_C592L3)

# Rename the columns for clarity
colnames(combined_counts_C592L) <- c("C591_1", "C591_2", "C591_3")

# Set the rownames to gene names or IDs if they are not already set
rownames(combined_counts_C592L) <- rownames(fcMUCK_C592L1$counts)  


##C593
fcMUCK_C593E1 <- featureCounts(c593e_bam_file1, annot.ext = annotation_file, isGTFAnnotationFile = TRUE,
                               GTF.featureType = "exon", GTF.attrType = "gene_id", 
                               useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)

fcMUCK_C593E2 <- featureCounts(c593e_bam_file2, annot.ext = annotation_file, isGTFAnnotationFile = TRUE,
                               GTF.featureType = "exon", GTF.attrType = "gene_id", 
                               useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)

fcMUCK_C593E3 <- featureCounts(c593e_bam_file3, annot.ext = annotation_file, isGTFAnnotationFile = TRUE,
                               GTF.featureType = "exon", GTF.attrType = "gene_id", 
                               useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)

# Extract the counts from each featureCounts output object
counts_C593E1 <- fcMUCK_C593E1$counts
counts_C593E2 <- fcMUCK_C593E2$counts
counts_C593E3 <- fcMUCK_C593E3$counts

# Combine the counts into a single dataframe
combined_counts_C593E <- cbind(counts_C593E1, counts_C593E2, counts_C593E3)

# Rename the columns for clarity
colnames(combined_counts_C593E) <- c("C59E_1", "C59E_2", "C59E_3")

# Set the rownames to gene names or IDs if they are not already set
rownames(combined_counts_C593E) <- rownames(fcMUCK_C593E1$counts)  

##C594
fcMUCK_C593L1 <- featureCounts(c593l_bam_file1, annot.ext = annotation_file, isGTFAnnotationFile = TRUE,
                               GTF.featureType = "exon", GTF.attrType = "gene_id", 
                               useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)

fcMUCK_C593L2 <- featureCounts(c593l_bam_file2, annot.ext = annotation_file, isGTFAnnotationFile = TRUE,
                               GTF.featureType = "exon", GTF.attrType = "gene_id", 
                               useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)

fcMUCK_C593L3 <- featureCounts(c593l_bam_file3, annot.ext = annotation_file, isGTFAnnotationFile = TRUE,
                               GTF.featureType = "exon", GTF.attrType = "gene_id", 
                               useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)

# Extract the counts from each featureCounts output object
counts_C593L1 <- fcMUCK_C593L1$counts
counts_C593L2 <- fcMUCK_C593L2$counts
counts_C593L3 <- fcMUCK_C593L3$counts

# Combine the counts into a single dataframe
combined_counts_C593L <- cbind(counts_C593L1, counts_C593L2, counts_C593L3)

# Rename the columns for clarity
colnames(combined_counts_C593L) <- c("C59L_1", "C59L_2", "C59L_3")

# Set the rownames to gene names or IDs if they are not already set
rownames(combined_counts_C593L) <- rownames(fcMUCK_C593L1$counts)  

#####

# Read the TSV files into dataframes with proper row names
combined_counts <- read.delim("combined_counts.tsv", sep = "\t", row.names = 1)
combined_counts_WNT <- read.delim("combined_countsWNT2.tsv", sep = "\t", row.names = 1)


# Convert row names to a column named "Gene" for each dataframe
combined_counts <- data.frame(Gene = rownames(combined_counts), combined_counts, row.names = NULL)
combined_counts_WNT <- data.frame(Gene = rownames(combined_counts_WNT), combined_counts_WNT, row.names = NULL)

combined_counts_C59E <- data.frame(Gene = rownames(combined_counts_C593E), combined_counts_C593E, row.names = NULL)
combined_counts_C59L <- data.frame(Gene = rownames(combined_counts_C593L), combined_counts_C593L, row.names = NULL)

# Merge dataframes for ES and WNT samples
combined_ES_WNT <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), 
                          list(combined_counts, combined_counts_WNT))

# Merge dataframes for ES and C59 samples
combined_ES_C59 <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), 
                          list(combined_counts, combined_counts_C59E, combined_counts_C59L))

# View the combined dataframes
print(combined_ES_WNT)
print(combined_ES_C59)

###################################################
# Convert the Gene column back to row names
row.names(combined_ES_WNT) <- combined_ES_WNT$Gene
combined_ES_WNT <- combined_ES_WNT[ , -which(names(combined_ES_WNT) == "Gene")]

row.names(combined_ES_C59) <- combined_ES_C59$Gene
combined_ES_C59 <- combined_ES_C59[ , -which(names(combined_ES_C59) == "Gene")]

# View the final dataframes
print(combined_ES_WNT)
print(combined_ES_C59)


# Remove rows (genes) with only zero values across all columns (samples)
cleaned_Wnt <- combined_ES_WNT[rowSums(combined_ES_WNT != 0) > 0, ]
cleaned_C59 <- combined_ES_C59[rowSums(combined_ES_C59 != 0) > 0, ]


# Optionally, if you want to save the cleaned matrix back to a TSV file
write.table(cleaned_Wnt, file = "cleaned_Wnt.tsv", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cleaned_C59, file = "cleaned_C59.tsv", sep = "\t", quote = FALSE, row.names = TRUE)


# Example col_data (replace with your actual metadata)
# Get sample names
sample_namesWNT <- colnames(cleaned_Wnt)
sample_namesC59 <-colnames(cleaned_C59)

# Define conditions
conditionsC59 <- factor(c(rep("control", 3), rep("treated", 6)))
conditionsWNT <- factor(c(rep("control", 3), rep("treated", 3)))

# Verify the conditions
conditionsC59
conditionsWNT


# Create col_data data frame
col_dataC59 <- data.frame(sampleName = sample_namesC59,
                          condition = conditionsC59,
                          stringsAsFactors = FALSE)
col_dataWNT <-data.frame(sampleName = sample_namesWNT,
                         condition = conditionsWNT,
                         stringsAsFactors = FALSE)

# Create DGEList object for data without zero value rows
dgeWNT <- DGEList(counts=cleaned_Wnt, group=conditionsWNT)
dgeC59 <- DGEList(counts=cleaned_C59, group=conditionsC59)


# Calculate normalized counts
norm_cleaned_countsWNT <- cpm(dgeWNT)
norm_cleaned_countsc59 <- cpm(dgeC59)

#Write to file
write.table(norm_cleaned_countsWNT, file="normalized_cleaned_WNT.tsv", sep="\t", row.names=TRUE, col.names=NA)
write.table(norm_cleaned_countsc59, file="normalized_cleaned_C59.tsv", sep="\t", row.names=TRUE, col.names=NA)

##########################################################

#####CEBPA and FOXA1 Datasets, Nagakawa et al. (2024)####
setwd("/Users/cassandradrazick/Desktop/BINF6999/SNagakawa")
data_dir <- "/Users/cassandradrazick/Desktop/BINF6999/SNagakawa/GSE244526_RAW"
# Load necessary libraries
library(readr)
library(dplyr)
library(tidyverse)
library(purrr)
library(tools)

# List all .genes.results.gz files
files <- list.files(data_dir, pattern = "*.genes.results.gz", full.names = TRUE)

# Check if files are found
if(length(files) == 0) {
  stop("No .genes.results.gz files found in the specified directory.")
}

# Initialize an empty list to store data frames
dfN <- list()

# Loop through each file and read data
for(file in files) {
  # Read the .gz file
  dfN <- read_tsv(file, col_names = TRUE)
  dfN <- as.data.frame(dfN)
  # Extract the sample name from the file name
  sample_name <- tools::file_path_sans_ext(basename(file))
  sample_name <- gsub(".genes.results", "", sample_name)
  
  # Select only the relevant columns 
  dfN <- dfN %>% select(gene_id, expected_count)
  
  # Rename the expected_count column to the sample name
  colnames(dfN)[2] <- sample_name
  
  # Add the data frame to the list
  dfN[[sample_name]] <- dfN
}

# Combine all data frames by gene_id
NagaDF <- purrr::reduce(dfN, full_join, by = "gene_id")

# Handle missing values by replacing NAs with 0
NagaDF[is.na(NagaDF)] <- 0
NagaDF <- as.data.frame(NagaDF)
# Set gene_id as row names
row.names(NagaDF) <- NagaDF$gene_id
NagaDF <- NagaDF[ , -1]  # Remove the gene_id column after setting row names


# Print the first few rows of the combined data frame to verify
print(head(NagaDF))

# Create a DGEList object for differential gene expression analysis
group <- factor(c(rep("Control", 3), rep("Treated,CEBPA,D1", 3), rep("Treated, CEBPA,D2", 3), rep("Treated, FOXA1,D1", 3), rep("Treated, FOXA1,D2", 3)))
dgeNAG <- DGEList(counts=NagaDF, group=group)

# Calculate normalized counts (CPM)
norm_countsNAG <- cpm(dgeNAG)

# View original column names to verify their format
print(colnames(norm_countsNAG))

# Use gsub to remove the unwanted part of the column names
colnames(norm_countsNAG) <- gsub("GSM\\d+_HuH-7_", "", colnames(norm_countsNAG))

# View modified column names to verify the changes
print(colnames(norm_countsNAG))

# Save normalized counts to a TSV file
write.table(norm_countsNAG, file=file.path(data_dir, "normalized_countsNAG.tsv"), sep="\t", row.names=TRUE, col.names=NA)

dfNAG <- read.table("/Users/cassandradrazick/Desktop/BINF6999/SNagakawa/normalized_countsNAG.tsv", header = TRUE, sep = "\t", row.names = 1)

# Check the structure of the dataframe
str(dfNAG)

#Remove first row
dfNAG <- dfNAG[-1]

##GO For this dataset - Gene names were included, lucky me :)
# Split row names
row_names <- rownames(dfNAG)
ensg_ids <- sub("\\..*$", "", row_names)  # Extract ENSG ID before the decimal
gene_names <- sub(".*_", "", row_names)    # Extract gene name after the underscore

# Create a new dataframe for ENSG IDs
dfnagENSG <- dfNAG
rownames(dfnagENSG) <- ensg_ids

# Create a new dataframe for gene names
dfnagGENES <- dfNAG
# Check for duplicate gene names
duplicated_gene_names <- gene_names[duplicated(gene_names)]
if(length(duplicated_gene_names) > 0) {
  message("Duplicate gene names found: ", paste(unique(duplicated_gene_names), collapse = ", "))
}

#Aggregate by taking the mean of duplicated rows
dfnagGENES_agg <- dfNAG %>%
  mutate(gene_name = gene_names) %>%
  group_by(gene_name) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

# Set row names and remove the gene_name column
dfnagGENES_agg<- as.data.frame(dfnagGENES_agg)
rownames(dfnagGENES_agg) <- dfnagGENES_agg$gene_name
dfnagGENES_agg <- dfnagGENES_agg[-1]

# Save the aggregated dataframe
write.table(dfnagGENES_agg, file = "dfnagGENES_agg.tsv", sep = "\t", row.names = TRUE, col.names = NA)
write.table(dfnagENSG, file = "dfnagENSG.tsv", sep = "\t", row.names = TRUE, col.names = NA)

# Read in the TSV files
dfnagGENES_agg <- read.table("dfnagGENES_agg.tsv", header = TRUE, sep = "\t", row.names = 1)

# Print original column names for reference
print(colnames(dfnagGENES_agg))

# Remove the prefix 'GSMXXXXXXX_HuH.7_' from column names
colnames(dfnagGENES_agg) <- gsub("^GSM[0-9]+_HuH\\.7_", "", colnames(dfnagGENES_agg))

# Print modified column names to verify
print(colnames(dfnagGENES_agg))

dfnagGENESC <- dfnagGENES_agg

# Print original column names for reference
print(colnames(dfnagGENES_agg))

# Create a vector of column names that should be kept for Control samples
control_columns <- grep("Control_[1-3]$", colnames(dfnagGENES_agg), value = TRUE)

# Create subsets for CEBPA and FOXA1 samples, keeping only Control columns

# Subset for CEBPA samples
cebp_cols <- grep("CEBPA", colnames(dfnagGENES_agg), value = TRUE)
cebp_subset <- dfnagGENES_agg[, c(control_columns, cebp_cols)]

# Subset for FOXA1 samples
foxa_cols <- grep("FOXA1", colnames(dfnagGENES_agg), value = TRUE)
foxa_subset <- dfnagGENES_agg[, c(control_columns, foxa_cols)]

# Print column names for the new subsets to verify
print(colnames(cebp_subset))
print(colnames(foxa_subset))

# Write cleaned data frames to TSV files
write.table(cebp_subset, "nagCEBPA.tsv", sep = "\t", col.names = NA, quote = FALSE)

write.table(foxa_subset, "nagFOXA1.tsv", sep = "\t", col.names = NA, quote = FALSE)

#Oops, forgot to remove 0s, now removing zeros
nagCEBPA <- read.table("nagCEBPA.tsv", header = TRUE, sep = "\t", row.names = 1 )
nagFOXA1 <- read.table("nagFOXA1.tsv", header = TRUE, sep = "\t", row.names = 1 )

# Remove rows containing only zeros from nagCEBPA
nagCEBPA_clean <- nagCEBPA[rowSums(nagCEBPA != 0) > 0, ]

# Remove rows containing only zeros from nagFOXA1
nagFOXA1_clean <- nagFOXA1[rowSums(nagFOXA1 != 0) > 0, ]

# Save the cleaned data frames to new TSV files
write.table(nagCEBPA_clean, file = "nagCEBPA_clean.tsv", sep = "\t", row.names = TRUE, col.names = NA)
write.table(nagFOXA1_clean, file = "nagFOXA1_clean.tsv", sep = "\t", row.names = TRUE, col.names = NA)

###################################################

#####CHIR99021 Dataset, Pagella et al. (2023)#####

library(BiocManager)
library(edgeR)
library(readr)
setwd("/Users/cassandradrazick/Desktop/BINF6999/CCantu")
# Read the raw counts matrix
count_data <- read.table("RNA_hESC_raw_counts_matrix.txt", header=TRUE, row.names=1, sep="\t")

# Replace NA values with zeros
count_data[is.na(count_data)] <- 0

# Define group factor (modify based on your experimental design)
group <- factor(c(rep("CHIR", 3), rep("CTRL", 3)))

# Identify rows (genes) with only zero values
non_zero_rows <- rowSums(count_data != 0) > 0

# Create a data frame without rows having only zero values
count_data_no_zeros <- count_data[non_zero_rows, ]

# Create DGEList object for data without zero value rows
dge_no_zeros <- DGEList(counts=count_data_no_zeros, group=group)
dge_no_zeros <- calcNormFactors(dge_no_zeros)

# Calculate normalized counts
norm_counts_no_zeros <- cpm(dge_no_zeros)

##########################################################
######Nkd1 Axin2 Dataset, Bell et al. (2024) #####
# Load packages
library(ShortRead)
library(Rsubread)
library(edgeR)
library(tidyverse)  
library(biomaRt)
setwd("/Users/cassandradrazick/Desktop/BINF6999/IanBell")
##Preprocess
##Feature counts 

# Define the annotation file 
annotation_file <- "/Users/cassandradrazick/Desktop/BINF6999/IanBell/Danio_rerio.GRCz11.112.gtf.gz"
# Define the .bam files
bam_files <- c(
  "double-mutant-50_S82_L001_Aligned.sortedByCoord.out.bam",
  "double-mutant-50_S82_L002_Aligned.sortedByCoord.out.bam",
  "double-mutant-50_S82_L003_Aligned.sortedByCoord.out.bam",
  "double-mutant-50_S82_L004_Aligned.sortedByCoord.out.bam",
  "Tu-50_S81_L001_Aligned.sortedByCoord.out.bam",
  "Tu-50_S81_L002_Aligned.sortedByCoord.out.bam",
  "Tu-50_S81_L003_Aligned.sortedByCoord.out.bam",
  "Tu-50_S81_L004_Aligned.sortedByCoord.out.bam")

# FeatureCounts to count reads
fcIAN <- featureCounts(bam_files, annot.ext = annotation_file, isGTFAnnotationFile = TRUE, 
                       GTF.featureType = "exon", GTF.attrType = "gene_id", 
                       useMetaFeatures = TRUE, isPairedEnd = TRUE)
count_matrixIAN <- fcIAN$counts

# Check the current column names
colnames(count_matrixIAN)

# Define the new column names in the specified order
new_column_names <- c("Wnt+ L001", "Wnt+ L002", "Wnt+ L003", "Wnt+ L004",
                      "Control L001", "Control L002", "Control L003", "Control L004")
# Assign the new column names to the count matrix
colnames(count_matrixIAN) <- new_column_names

# Verify the changes
colnames(count_matrixIAN)

# Remove rows (genes) with only zero values across all columns (samples)
cleaned_countsIAN <- count_matrixIAN[rowSums(count_matrixIAN != 0) > 0, ]


# Save the cleaned matrix back to a TSV file
write.table(cleaned_countsIAN, file = "cleaned_count_matrixIAN.tsv", sep = "\t", quote = FALSE, row.names = TRUE)
setwd("/Users/cassandradrazick/Desktop/BINF6999/IanBell/")

# Get sample names
sample_names <- colnames(count_matrixIAN)

# Define conditions: first 4 samples treated, last 4 samples control
conditions <- c(rep("treated", 4), rep("control", 4))

# Verify the conditions
conditions
# Assuming all samples are treated

# Create col_data data frame
col_data <- data.frame(sampleName = sample_names,
                       condition = conditions,
                       stringsAsFactors = FALSE)


groupIAN <- factor(c(rep("Wnt+", 4), rep("WT", 4)))

# Create DGEList object for data without zero value rows
dgeIAN <- DGEList(counts=cleaned_countsIAN, group=groupIAN)
dgeIAN <- calcNormFactors(dgeIAN)

# Calculate normalized counts
norm_cleaned_countsIAN <- cpm(dgeIAN)

#Write to file
write.table(norm_cleaned_countsIAN, file="normalized_cleaned_IAN.tsv", sep="\t", row.names=TRUE, col.names=NA)

counts_matrixIAN <- read.table("normalized_cleaned_IAN.tsv")





