import argparse
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr, linregress
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np

def main():
    parser = argparse.ArgumentParser(description='Scatter plot with regression line for two datasets.')
    parser.add_argument('-i', '--input1', required=True, help='Path to the first CSV file')
    parser.add_argument('-j', '--input2', required=True, help='Path to the second CSV file')
    parser.add_argument('-x', '--xlabel', required=True, help='Label for the X-axis')
    parser.add_argument('-y', '--ylabel', required=True, help='Label for the Y-axis')
    parser.add_argument('-o', '--output', required=True, help='Path to save the output image file')
    args = parser.parse_args()

    # Read the CSV files
    df1 = pd.read_csv(args.input1, header=None, names=['pair', 'value1'])
    df2 = pd.read_csv(args.input2, header=None, names=['pair', 'value2'])

    # Merge the two dataframes on 'pair'
    merged_df = pd.merge(df1, df2, on='pair')

    # Extract the values
    values1 = merged_df['value1']
    values2 = merged_df['value2']

    # Pearson correlation coefficient
    pearson_corr, _ = pearsonr(values1, values2)
    print(f"Pearson correlation coefficient: {pearson_corr}")

    # Spearman correlation coefficient
    spearman_corr, _ = spearmanr(values1, values2)
    print(f"Spearman correlation coefficient: {spearman_corr}")

    # Perform linear regression
    slope, intercept, r_value, p_value, std_err = linregress(values1, values2)
    r_squared = r_value ** 2
    print(f"R-squared: {r_squared}")

    # Mean Squared Error (MSE)
    mse = mean_squared_error(values1, values2)
    print(f"Mean Squared Error (MSE): {mse}")

    # Mean Absolute Error (MAE)
    mae = mean_absolute_error(values1, values2)
    print(f"Mean Absolute Error (MAE): {mae}")

    # Cosine Similarity
    cosine_sim = cosine_similarity(values1.values.reshape(1, -1), values2.values.reshape(1, -1))
    print(f"Cosine Similarity: {cosine_sim[0][0]}")

    # Scatter Plot with Regression Line
    plt.figure(figsize=(8, 6))
    plt.scatter(values1, values2, alpha=0.02, marker='+', color='blue')  # Use small dots
    
    # Add the regression line
    plt.plot(values1, slope * values1 + intercept, color='red', label=f'Regression Line (R²={r_squared:.3f})')
    
    # Add x=y reference line
    plt.plot([min(values1), max(values1)], [min(values1), max(values1)], color='black', linestyle='--', label='x = y')

    # Add Pearson and Spearman correlation to the legend
    plt.legend([
        'x = y',
        f'Regression Line (R²={r_squared:.3f})',
        f'Pearson: {pearson_corr:.3f}',
        f'Spearman: {spearman_corr:.3f}'
    ])
    
    plt.title('Scatter Plot with Regression Line')
    plt.xlabel(args.xlabel)
    plt.ylabel(args.ylabel)
    
    # Save scatter plot
    plt.savefig(args.output, format='png')
    plt.show()

if __name__ == "__main__":
    main()

