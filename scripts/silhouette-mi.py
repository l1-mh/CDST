import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score, adjusted_mutual_info_score
from sklearn.cluster import AgglomerativeClustering
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import argparse
import os
from joblib import Parallel, delayed
import logging
import time

# === Configure logging system ===
logging.basicConfig(
    format="%(asctime)s - %(message)s",
    level=logging.INFO,
    datefmt="%H:%M:%S"
)

# === MI computation task (parallelizable) ===
def compute_mi(i, j, labels_list):
    if i == j:
        return i, j, 1.0
    else:
        mi = adjusted_mutual_info_score(labels_list[i], labels_list[j])
        return i, j, mi


def analyze_clustering(input_file, output_dir, step_size, mode, threads):
    # === Load distance matrix ===
    distance_matrix = pd.read_csv(input_file, index_col=0)

    # Make matrix symmetric by keeping the minimum value
    distance_matrix = np.minimum(distance_matrix, distance_matrix.T)

    # Set diagonal to 0
    np.fill_diagonal(distance_matrix.values, 0)

    os.makedirs(output_dir, exist_ok=True)

    labels_list = []
    thresholds = []
    silhouette_scores = []

    # === Perform hierarchical clustering by distance or by number of clusters ===
    if mode == 'dist':
        max_distance = np.max(distance_matrix.values)
        thresholds = np.arange(0, max_distance + step_size, step_size)

        logging.info(f"Max distance: {max_distance:.3f}, Step size: {step_size:.3f}")

        for threshold in thresholds:
            clustering = AgglomerativeClustering(
                n_clusters=None,
                metric='precomputed',
                linkage='average',
                distance_threshold=threshold
            )
            labels = clustering.fit_predict(distance_matrix)
            labels_list.append(labels)

            num_clusters = len(set(labels))
            if num_clusters == 1:
                logging.info(f"All samples merged at threshold = {threshold:.3f}")
                break

            if 2 <= num_clusters <= len(distance_matrix) - 1:
                score = silhouette_score(distance_matrix, labels, metric='precomputed')
                silhouette_scores.append([threshold, score])
                logging.info(f"Threshold: {threshold:.3f}, Silhouette Score: {score:.3f}")
            else:
                silhouette_scores.append([threshold, None])

    elif mode == 'nclust':
        min_clusters = 2
        max_clusters = len(distance_matrix)
        thresholds = list(range(min_clusters, max_clusters + 1, int(step_size)))

        logging.info(f"Clustering from {min_clusters} to {max_clusters} clusters with step size = {step_size}")

        for n_clusters in thresholds:
            clustering = AgglomerativeClustering(
                n_clusters=n_clusters,
                metric='precomputed',
                linkage='average'
            )
            labels = clustering.fit_predict(distance_matrix)
            labels_list.append(labels)

            num_clusters = len(set(labels))
            if num_clusters == 1:
                logging.info(f"All samples merged at {n_clusters} clusters.")
                break

            if 2 <= num_clusters <= len(distance_matrix) - 1:
                score = silhouette_score(distance_matrix, labels, metric='precomputed')
                silhouette_scores.append([n_clusters, score])
                logging.info(f"Clusters: {n_clusters}, Silhouette Score: {score:.3f}")
            else:
                silhouette_scores.append([n_clusters, None])
    else:
        raise ValueError("Invalid mode: Use 'dist' or 'nclust'")

    # Adjust threshold list to match number of clustering results
    thresholds = thresholds[:len(labels_list)]

    # === Save silhouette score results ===
    silhouette_df = pd.DataFrame(silhouette_scores, columns=["Threshold/Clusters", "Silhouette Score"])
    silhouette_file = os.path.join(output_dir, "silhouette_scores.csv")
    silhouette_df.to_csv(silhouette_file, index=False)
    logging.info(f"Silhouette scores saved to: {silhouette_file}")

    # === Plot silhouette score line chart ===
    plt.figure(figsize=(8, 6))
    plt.plot(
        silhouette_df["Threshold/Clusters"],
        silhouette_df["Silhouette Score"],
        linestyle='-',  # Line only, no data markers
        color='blue'
    )
    plt.title(f"Silhouette Score vs. {'Distance Threshold' if mode == 'dist' else 'Number of Clusters'}")
    plt.xlabel("Threshold/Clusters")
    plt.ylabel("Silhouette Score")
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, "silhouette_plot.png"))
    logging.info(f"Silhouette plot saved to: silhouette_plot.png")

    # === Compute mutual information matrix ===
    n = len(labels_list)
    logging.info("Starting MI calculation...")
    start_time = time.time()

    results = Parallel(n_jobs=threads, batch_size=1)(
        delayed(compute_mi)(i, j, labels_list)
        for i in range(n) for j in range(i, n)
    )

    elapsed_time = time.time() - start_time
    logging.info(f"MI calculation completed in {elapsed_time:.2f} seconds.")

    # === Build MI matrix from results ===
    mi_matrix = np.zeros((n, n))
    for i, j, mi in results:
        mi_matrix[i][j] = mi_matrix[j][i] = mi

    mi_matrix_df = pd.DataFrame(mi_matrix, index=thresholds, columns=thresholds)
    mi_matrix_file = os.path.join(output_dir, "mi_matrix.csv")
    mi_matrix_df.to_csv(mi_matrix_file)
    logging.info(f"Mutual Information matrix saved to: {mi_matrix_file}")

    # === Plot MI heatmap with color midpoint at 0.9 and ticks at 0.1 intervals ===
    plt.figure(figsize=(10, 8))
    norm = TwoSlopeNorm(vmin=0, vcenter=0.9, vmax=1)

    im = plt.imshow(mi_matrix, cmap='coolwarm', interpolation='nearest', norm=norm)
    cbar = plt.colorbar(im)
    cbar.set_ticks(np.arange(0, 1.1, 0.1))
    cbar.set_label('Adjusted Mutual Information')

    plt.title(f'Mutual Information Matrix ({mode})')

    # Show only 5 evenly spaced ticks on axes
    num_ticks = 5
    x_ticks = np.linspace(0, len(thresholds) - 1, num_ticks).astype(int)
    x_labels = [f'{thresholds[i]:.2f}' for i in x_ticks]

    plt.xticks(ticks=x_ticks, labels=x_labels, rotation=45)
    plt.yticks(ticks=x_ticks, labels=x_labels)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "mi_matrix.png"))
    logging.info(f"Mutual Information plot saved to: mi_matrix.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze clustering consistency using Silhouette and Mutual Information.")
    parser.add_argument('-i', '--input', type=str, required=True, help="Input distance matrix file (CSV)")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output directory for results")
    parser.add_argument('-s', '--step', type=float, default=1.0, help="Step size for thresholds or cluster counts")
    parser.add_argument('-m', '--mode', type=str, choices=['dist', 'nclust'], required=True, help="Clustering mode: 'dist' or 'nclust'")
    parser.add_argument('-t', '--threads', type=int, default=-1, help="Number of parallel threads to use")

    args = parser.parse_args()
    analyze_clustering(args.input, args.output, args.step, args.mode, args.threads)
