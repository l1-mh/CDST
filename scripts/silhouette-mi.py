from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score, adjusted_mutual_info_score
from sklearn.cluster import AgglomerativeClustering
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from joblib import Parallel, delayed
import argparse
import logging
import time


# --------------------------- Logging Configuration ---------------------------
logging.basicConfig(
    format="%(asctime)s - %(message)s",
    level=logging.INFO,
    datefmt="%H:%M:%S"
)

# --------------------------- Utility Functions ---------------------------

def load_distance_matrix(file_path):
    matrix = pd.read_csv(file_path, index_col=0)
    matrix = np.minimum(matrix, matrix.T)
    np.fill_diagonal(matrix.values, 0)
    return matrix


def compute_mi(i, j, labels_list):
    if i == j:
        return i, j, 1.0
    mi = adjusted_mutual_info_score(labels_list[i], labels_list[j])
    return i, j, mi


def save_dataframe(df, output_path, name):
    path = output_path / name
    df.to_csv(path, index=False)
    logging.info(f"{name} saved to: {path}")


def run_clustering(matrix, mode, step, start, end):
    labels_list = []
    scores = []

    if mode == 'dist':
        max_dist = np.max(matrix.values)
        start = 0 if start is None else start
        end = max_dist if end is None else end
        thresholds = np.arange(start, end + step, step)

        logging.info(f"Running clustering by distance from {start} to {end} (step {step})")

        for threshold in thresholds:
            model = AgglomerativeClustering(n_clusters=None, metric='precomputed',
                                            linkage='average', distance_threshold=threshold)
            labels = model.fit_predict(matrix)
            labels_list.append(labels)

            n_clusters = len(set(labels))
            if n_clusters == 1:
                logging.info(f"All samples merged at threshold = {threshold:.3f}")
                break

            score = silhouette_score(matrix, labels, metric='precomputed') if 2 <= n_clusters <= len(matrix) - 1 else None
            scores.append([threshold, score, n_clusters])
            logging.info(f"Threshold: {threshold:.3f}, Silhouette Score: {score}, Clusters: {n_clusters}")

    else:
        max_clusters = len(matrix)
        start = 2 if start is None else int(start)
        end = max_clusters if end is None else int(end)
        thresholds = list(range(start, end + 1, int(step)))

        logging.info(f"Running clustering by number of clusters from {start} to {end} (step {step})")

        for n_clusters in thresholds:
            model = AgglomerativeClustering(n_clusters=n_clusters, metric='precomputed', linkage='average')
            labels = model.fit_predict(matrix)
            labels_list.append(labels)

            score = silhouette_score(matrix, labels, metric='precomputed') if 2 <= n_clusters <= len(matrix) - 1 else None
            scores.append([n_clusters, score, np.nan])
            logging.info(f"Clusters: {n_clusters}, Silhouette Score: {score}")

    return labels_list, scores, thresholds[:len(labels_list)]


def plot_silhouette(scores_df, mode, output_path, img_format):
    x = scores_df.iloc[:, 0]
    sil = scores_df["Silhouette Score"]
    aux = scores_df.iloc[:, 2]

    fig, ax1 = plt.subplots(figsize=(8, 6))
    ax1.plot(x, sil, color='blue')
    ax1.set_xlabel(scores_df.columns[0])
    ax1.set_ylabel("Silhouette Score", color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')

    ax2 = ax1.twinx()
    ax2.plot(x, aux, linestyle='--', color='red')
    ax2.set_ylabel(scores_df.columns[2], color='red')
    ax2.tick_params(axis='y', labelcolor='red')

    plt.title("Silhouette Score and " + scores_df.columns[2])
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path / f"silhouette_plot.{img_format}", format=img_format)
    logging.info(f"Silhouette plot saved as {img_format}.")


def compute_mi_matrix(labels_list, threads, thresholds, output_path, img_format):
    n = len(labels_list)
    logging.info("Starting Mutual Information matrix calculation...")

    start = time.time()
    results = Parallel(n_jobs=threads, batch_size=1)(
        delayed(compute_mi)(i, j, labels_list)
        for i in range(n) for j in range(i, n)
    )
    elapsed = time.time() - start
    logging.info(f"MI calculation completed in {elapsed:.2f} seconds.")

    mi_matrix = np.zeros((n, n))
    for i, j, val in results:
        mi_matrix[i][j] = mi_matrix[j][i] = val

    df = pd.DataFrame(mi_matrix, index=thresholds, columns=thresholds)
    df.to_csv(output_path / "mi_matrix.csv")
    logging.info("MI matrix saved.")

    plt.figure(figsize=(10, 8))
    norm = TwoSlopeNorm(vmin=0, vcenter=0.9, vmax=1)
    im = plt.imshow(mi_matrix, cmap='coolwarm', interpolation='nearest', norm=norm)
    cbar = plt.colorbar(im)
    cbar.set_ticks(np.arange(0, 1.1, 0.1))
    cbar.set_label("Adjusted Mutual Information")

    ticks = np.linspace(0, len(thresholds) - 1, 5).astype(int)
    tick_labels = [f'{thresholds[i]:.2f}' for i in ticks]
    plt.xticks(ticks, tick_labels, rotation=45)
    plt.yticks(ticks, tick_labels)
    plt.title("Mutual Information Matrix")
    plt.tight_layout()
    plt.savefig(output_path / f"mi_matrix.{img_format}", format=img_format)
    logging.info(f"Mutual Information heatmap saved as {img_format}.")


# --------------------------- Main Function ---------------------------

def main(args):
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)

    matrix = load_distance_matrix(args.input)
    labels_list, scores, thresholds = run_clustering(matrix, args.mode, args.step, args.start, args.end)

    col_names = ["Threshold" if args.mode == 'dist' else "Clusters",
                 "Silhouette Score",
                 "Clusters" if args.mode == 'dist' else "Distance"]

    scores_df = pd.DataFrame(scores, columns=col_names)
    save_dataframe(scores_df, output_path, "silhouette_scores.csv")
    plot_silhouette(scores_df, args.mode, output_path, args.format)

    if args.mi:
        compute_mi_matrix(labels_list, args.threads, thresholds, output_path, args.format)


# --------------------------- CLI ---------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clustering evaluation using silhouette and mutual information.")
    parser.add_argument("-i", "--input", required=True, help="Distance matrix (CSV)")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("-s", "--step", type=float, default=1.0, help="Step size")
    parser.add_argument("-m", "--mode", choices=["dist", "nclust"], required=True, help="Clustering mode")
    parser.add_argument("-t", "--threads", type=int, default=-1, help="Threads for MI computation")
    parser.add_argument("--start", type=float, help="Start value")
    parser.add_argument("--end", type=float, help="End value")
    parser.add_argument("-M", "--mi", action="store_true", help="If set, compute Mutual Information matrix")
    parser.add_argument("--format", type=str, default="jpg", choices=["jpg", "png", "tif", "pdf", "svg"],
                        help="Image format for output plots (default: jpg)")

    args = parser.parse_args()
    main(args)

