import pandas as pd
import numpy as np
import argparse
import os
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
import matplotlib.pyplot as plt

def evaluate_clustering(matrix_file, label_file, mode, step, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # === Load distance matrix ===
    dist_df = pd.read_csv(matrix_file, index_col=0)

    # === Load and clean reference labels ===
    ref_df = pd.read_csv(label_file)

    # Drop rows with missing Cluster labels
    ref_df = ref_df.dropna(subset=["Cluster"])

    # Ensure sample names exist in the matrix
    valid_samples = list(set(ref_df["Sample"]) & set(dist_df.index))
    ref_df = ref_df.set_index("Sample").loc[valid_samples]
    ref_labels = ref_df["Cluster"].values

    # Filter distance matrix to valid samples
    dist_df = dist_df.loc[valid_samples, valid_samples]
    sample_names = dist_df.index.tolist()
    dist_matrix = np.minimum(dist_df.values, dist_df.values.T)
    np.fill_diagonal(dist_matrix, 0)

    # === Clustering and score calculation ===
    results_ari = []
    results_ami = []

    if mode == 'dist':
        thresholds = np.arange(0, np.max(dist_matrix), step)
        for t in thresholds:
            model = AgglomerativeClustering(n_clusters=None, distance_threshold=t, linkage='average', metric='precomputed')
            pred_labels = model.fit_predict(dist_matrix)
            ari = adjusted_rand_score(ref_labels, pred_labels)
            ami = adjusted_mutual_info_score(ref_labels, pred_labels)
            results_ari.append((t, ari))
            results_ami.append((t, ami))
            print(f"Threshold={t:.2f}: ARI={ari:.3f}, AMI={ami:.3f}")
    else:  # nclust
        thresholds = range(2, len(ref_labels))
        for n in thresholds:
            model = AgglomerativeClustering(n_clusters=n, linkage='average', metric='precomputed')
            pred_labels = model.fit_predict(dist_matrix)
            ari = adjusted_rand_score(ref_labels, pred_labels)
            ami = adjusted_mutual_info_score(ref_labels, pred_labels)
            results_ari.append((n, ari))
            results_ami.append((n, ami))
            print(f"Clusters={n}: ARI={ari:.3f}, AMI={ami:.3f}")

    # === Save results ===
    ari_df = pd.DataFrame(results_ari, columns=["Threshold/Clusters", "ARI"])
    ami_df = pd.DataFrame(results_ami, columns=["Threshold/Clusters", "AMI"])
    ari_df.to_csv(os.path.join(output_dir, "ari_scores.csv"), index=False)
    ami_df.to_csv(os.path.join(output_dir, "ami_scores.csv"), index=False)

    # === Plot ===
    x_vals, y_ari = zip(*results_ari)
    _, y_ami = zip(*results_ami)

    plt.figure(figsize=(8, 6))
    plt.plot(x_vals, y_ari, label="Adjusted Rand Index", color="darkgreen", linestyle="-")
    plt.plot(x_vals, y_ami, label="Adjusted Mutual Information", color="darkorange", linestyle="-")
    plt.xlabel("Distance Threshold" if mode == 'dist' else "Number of Clusters")
    plt.ylabel("Score")
    plt.title("Clustering Consistency with Reference Labels")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plot_path = os.path.join(output_dir, "ari_ami_plot.png")
    plt.savefig(plot_path, dpi=300)
    print(f"\nPlot saved to: {plot_path}")

# === CLI entry ===
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate clustering consistency using ARI and AMI against reference labels.")
    parser.add_argument('-i', '--input', required=True, help="Input distance matrix (CSV)")
    parser.add_argument('-l', '--labels', required=True, help="Reference label file (CSV with 'Sample' and 'Cluster')")
    parser.add_argument('-m', '--mode', choices=['dist', 'nclust'], required=True, help="Clustering mode: 'dist' or 'nclust'")
    parser.add_argument('-s', '--step', type=float, default=50.0, help="Step size for threshold or number of clusters (default: 50)")
    parser.add_argument('-o', '--output', required=True, help="Output directory to save results")

    args = parser.parse_args()
    evaluate_clustering(args.input, args.labels, args.mode, args.step, args.output)
