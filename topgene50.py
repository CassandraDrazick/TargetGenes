import pandas as pd
import numpy as np
import sys
from scipy.stats import ttest_ind, zscore
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import LeaveOneOut
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
import matplotlib.pyplot as plt

def main(expression_file, conditions_file, output_file_prefix):
    # Load the gene expression data from a tab-separated file
    data = pd.read_csv(expression_file, sep='\t', index_col=0)
    
    # Load the conditions file which contains sample labels
    conditions = pd.read_csv(conditions_file, sep='\t')
    
    # Identify treated and control samples based on the conditions file
    treated_samples = conditions[conditions['Label'] == 'Treated']['Sample'].tolist()
    control_samples = conditions[conditions['Label'] == 'Control']['Sample'].tolist()
    
    # Combine all samples into one list and create a corresponding list of labels
    all_samples = treated_samples + control_samples
    labels = ['Treated'] * len(treated_samples) + ['Control'] * len(control_samples)
    
    # Ensure the data columns match the sample names
    data = data[all_samples]  # Data is transposed to have genes as rows
    X = data.T.values  # Feature matrix
    y = np.array(labels)  # Labels for classification
    
    # Perform PCA on the features
    pca = PCA(n_components=min(X.shape) - 1)
    pca.fit(X)
    
    # Retrieve PCA loadings and create a DataFrame of loadings
    loadings = pca.components_
    loadings_df = pd.DataFrame(loadings.T, index=data.index, columns=[f'PC{i+1}' for i in range(loadings.shape[0])])
    
    # Calculate the absolute value of loadings for each gene and sum them to determine importance
    abs_loadings_df = loadings_df.abs()
    importance = abs_loadings_df.sum(axis=1)
    
    # Get top features based on overall importance
    top_genes_pca = importance.sort_values(ascending=False).head(50)
    top_genes_pca_list = top_genes_pca.index.tolist()
    top_genes_pca.to_csv(f"{output_file_prefix}_top_features_from_pca.csv", header=True)
    
    # Plot top features based on the first principal component (PC1)
    top_genes_pc1 = loadings_df.loc[top_genes_pca_list, 'PC1'].abs().sort_values(ascending=False)
    plt.figure(figsize=(12, 8))
    top_genes_pc1.plot(kind='bar', color='skyblue')
    plt.xlabel('Gene')
    plt.ylabel('Absolute Loading')
    plt.title(f'Top Genes Based on PC1 Loadings for {output_file_prefix} Dataset')
    plt.xticks(rotation=90)
    plt.savefig(f"{output_file_prefix}_top_genes_pc1.png")
    plt.show()
    
    # Convert labels to numeric values for classification
    le = LabelEncoder()
    y_numeric = le.fit_transform(y)
    
    # Initialize Leave-One-Out Cross-Validation (LOOCV) and machine learning models
    loo = LeaveOneOut()
    models = {
        "SVM": SVC(kernel='linear', probability=True, random_state=444),
        "RandomForest": RandomForestClassifier(random_state=444),
        "NeuralNetwork": MLPClassifier(random_state=444, max_iter=1000)
    }
    
    classification_results = []  # To store classification metrics
    top_genes_dict = {"SVM": [], "RandomForest": [], "NeuralNetwork": []}  # To store top genes from each model
    
    for train_index, test_index in loo.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y_numeric[train_index], y_numeric[test_index]

        for model_name, model in models.items():
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)
            y_pred_proba = model.predict_proba(X_test)[:, 1] if hasattr(model, "predict_proba") else y_pred

            # Classification metrics
            acc = accuracy_score(y_test, y_pred)
            prec = precision_score(y_test, y_pred, average='weighted', zero_division=0)
            rec = recall_score(y_test, y_pred, average='weighted', zero_division=0)
            f1 = f1_score(y_test, y_pred, average='weighted', zero_division=0)
            classification_results.append({
                "Model": model_name,
                "Accuracy": acc,
                "Precision": prec,
                "Recall": rec,
                "F1 Score": f1
            })

            # Extract feature importances from the models
            if model_name == "RandomForest" and hasattr(model, 'feature_importances_'):
                feature_importances = model.feature_importances_
            elif model_name == "SVM" and hasattr(model, 'coef_'):
                feature_importances = model.coef_[0]
            elif model_name == "NeuralNetwork" and hasattr(model, 'coefs_'):
                feature_importances = np.sum(np.abs(model.coefs_[0]), axis=1)
            else:
                feature_importances = np.zeros(X_train.shape[1])

            # Get top 50 genes based on feature importance
            top_indices = np.argsort(np.abs(feature_importances))[::-1][:50]
            top_genes_list = data.index[top_indices].tolist()
            top_genes_dict[model_name].extend(top_genes_list)
    
    # Save machine learning model results to a CSV file
    combined_top_genes_df = pd.DataFrame({model: pd.Series(genes).head(50).values for model, genes in top_genes_dict.items()})
    combined_top_genes_df.index = range(1, len(combined_top_genes_df) + 1)
    
    # Add PCA top genes to the combined DataFrame
    combined_top_genes_df['PCA'] = pd.Series(top_genes_pca_list).head(50).values
    combined_top_genes_df.to_csv(f"{output_file_prefix}_top50_genes_combined.csv", index_label='Rank')

    # Save mean classification results
    classification_results_df = pd.DataFrame(classification_results)
    average_metrics = classification_results_df.groupby('Model').mean(numeric_only=True).reset_index()
    average_metrics.to_csv(f"{output_file_prefix}_average_classification_metrics.csv", index=False)
    
    # Perform t-tests for each gene to compare treated vs. control
    gene_names = []
    p_values = []
    t_values = []
    effect_sizes = [] 
    for gene in data.index:
        treated = data.loc[gene, treated_samples].dropna().values
        control = data.loc[gene, control_samples].dropna().values

        if len(treated) < 2 or len(control) < 2 or treated.std() == 0 or control.std() == 0:
            continue

        t_stat, p_val = ttest_ind(treated, control, equal_var=False)
        pooled_std = ((treated.std()**2 + control.std()**2) / 2)**0.5
        effect_size = (treated.mean() - control.mean()) / pooled_std

        gene_names.append(gene)
        p_values.append(p_val)
        t_values.append(t_stat)
        effect_sizes.append(effect_size)

    # Correct p-values for multiple testing
    _, p_values_corrected, _, _ = multipletests(p_values, method='fdr_bh')

    # Create a DataFrame for t-test results
    ttest_results_df = pd.DataFrame({
        'Gene': gene_names,
        'p-value': p_values,
        'corrected p-value': p_values_corrected,
        't-statistic': t_values,
        'effect size': effect_sizes
    })
    
    # Compute z-scores for t-statistics
    ttest_results_df['z-score'] = zscore(ttest_results_df['t-statistic'])

    # Identify top 50 genes based on corrected p-value
    top_genes_pval = ttest_results_df.nsmallest(50, 'corrected p-value')

    # Sort these top 50 genes by z-score
    top_genes_ttest = top_genes_pval.sort_values(by='z-score', ascending=False)

    # Save t-test results for the top 50 genes
    top_genes_ttest.to_csv(f"{output_file_prefix}_ttests_results.csv", sep='\t', index=False)
    
    # Combine top genes from PCA, machine learning models, and t-tests into one DataFrame
    combined_top_genes_list = {
        "Top 50 from PCA": top_genes_pca_list,
        "Top 50 from SVM": pd.Series(top_genes_dict["SVM"]).head(50).tolist(),
        "Top 50 from RandomForest": pd.Series(top_genes_dict["RandomForest"]).head(50).tolist(),
        "Top 50 from NeuralNetwork": pd.Series(top_genes_dict["NeuralNetwork"]).head(50).tolist(),
        "Top 50 from t-tests": top_genes_ttest['Gene'].tolist()
    }
    
    combined_top_genes_all_df = pd.DataFrame(combined_top_genes_list)
    combined_top_genes_all_df.index = range(1, 51)
    combined_top_genes_all_df.to_csv(f"{output_file_prefix}_combined_top_genes_all_sources.csv", index_label='Rank')

if __name__ == "__main__":
    # Check if the correct number of command-line arguments are provided
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <expression_file> <conditions_file> <output_file_prefix>")
        sys.exit(1)
    
    # Get command-line arguments
    expression_file = sys.argv[1]
    conditions_file = sys.argv[2]
    output_file_prefix = sys.argv[3]
    
    # Run the main function with provided arguments
    main(expression_file, conditions_file, output_file_prefix)
