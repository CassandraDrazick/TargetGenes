import pandas as pd
import sys
import os
from collections import defaultdict
#Define a function to merge gene counts into a single matrix
def create_gene_count_matrix(input_files, output_file_prefix):
    gene_counts = defaultdict(lambda: defaultdict(int))
    gene_occurrences = defaultdict(int)
    
    # Determine the model columns for each file and maintain the order
    file_model_columns = {}
    
    for file in input_files:
        dataset_name = os.path.basename(file).split('_')[0]
        df = pd.read_csv(file, sep=',', index_col=0)
        columns = [f"{dataset_name} {col.split(' ')[-1]}" for col in df.columns]
        file_model_columns[file] = columns
    
    # Process each file to count occurrences of each gene across models
    for file in input_files:
        dataset_name = os.path.basename(file).split('_')[0]
        df = pd.read_csv(file, sep=',', index_col=0)

        for column in df.columns:
            model = column.split(' ')[-1]  # Extract model name from column
            column_name = f"{dataset_name} {model}"
            if column_name not in file_model_columns[file]:
                continue

            genes = df[column].dropna().tolist()  # Drop NaNs and get gene list 
            for rank, gene in enumerate(genes, start=1):
                if gene:  # Avoid counting empty strings or NaNs
                    gene_counts[gene][column_name] += 1
                    gene_occurrences[gene] += 1
    
    # Create a DataFrame for the gene counts matrix
    columns = ['Gene', 'Total Occurrences']
    
    # Add model columns in the order of input files
    for file in input_files:
        columns.extend(file_model_columns[file])
    
    data = []
    
    for gene, counts in gene_counts.items():
        if gene_occurrences[gene] >= 2:  # Include genes that appear 2 or more times
            row = [gene, gene_occurrences[gene]]
            for col in columns[2:]:
                count = counts.get(col, 0)
                row.append(count)
            data.append(row)
    
    # Create DataFrame
    df = pd.DataFrame(data, columns=columns)
    
    # Sort DataFrame by 'Total Occurrences'
    df = df.sort_values(by='Total Occurrences', ascending=False)
    
    # Save to CSV
    df.to_csv(f"{output_file_prefix}_gene_count_matrix.csv", index=False)
    print(f"Gene count matrix saved to {output_file_prefix}_gene_count_matrix.csv")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <output_file_prefix> <input_file1> [<input_file2> ...]")
        sys.exit(1)
    #Usage of sys.argv to use with other datasets
    output_file_prefix = sys.argv[1]
    input_files = sys.argv[2:]
    create_gene_count_matrix(input_files, output_file_prefix)
