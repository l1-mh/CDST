# compare_distance_lists.py
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import pearsonr, spearmanr, linregress, t
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.metrics.pairwise import cosine_similarity
from pathlib import Path
import numpy as np
import sys

def infer_output_path(output_arg: str, img_format: str) -> str:
    p = Path(output_arg)
    if p.suffix.lower() != f".{img_format}":
        p = p.with_suffix(f".{img_format}")
    p.parent.mkdir(parents=True, exist_ok=True)
    return str(p)

def log(msg: str):
    sys.stdout.write(msg + "\n")
    sys.stdout.flush()

def auto_metrics_location(values1, values2):
    """自动选择散点稀疏的象限"""
    x_med = np.median(values1)
    y_med = np.median(values2)
    quadrants = {
        "upper right": 0,
        "upper left": 0,
        "lower left": 0,
        "lower right": 0
    }
    for x, y in zip(values1, values2):
        if x >= x_med and y >= y_med:
            quadrants["upper right"] += 1
        elif x < x_med and y >= y_med:
            quadrants["upper left"] += 1
        elif x < x_med and y < y_med:
            quadrants["lower left"] += 1
        else:
            quadrants["lower right"] += 1
    return min(quadrants, key=quadrants.get)

def main():
    parser = argparse.ArgumentParser(
        description="Scatter plot with regression line and metrics box for two datasets."
    )
    parser.add_argument("-i", "--input1", required=True)
    parser.add_argument("-j", "--input2", required=True)
    parser.add_argument("-x", "--xlabel", required=True)
    parser.add_argument("-y", "--ylabel", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--format", default="pdf",
                        choices=["pdf", "png", "tif", "tiff", "svg", "jpg", "jpeg"])
    parser.add_argument("--stats-loc", default=None,
                        choices=["upper left","upper right","lower left","lower right","center left","center right","center", None],
                        help="Metrics legend location; if not set, auto-choose least crowded quadrant")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--marker", default="+")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--ci", type=float, help="Confidence interval percentage, e.g., 95 for 95% CI (default: None)")
    args = parser.parse_args()

    # --- IO ---
    log(f"Reading input files: {args.input1} , {args.input2}")
    df1 = pd.read_csv(args.input1, header=0, names=["pair", "value1"])
    df2 = pd.read_csv(args.input2, header=0, names=["pair", "value2"])

    log("Merging on column 'pair' ...")
    merged_df = pd.merge(df1, df2, on="pair")
    n_pairs = merged_df.shape[0]
    log(f"Merged pairs: {n_pairs}")

    values1 = merged_df["value1"].astype(float)
    values2 = merged_df["value2"].astype(float)

    # --- Metrics ---
    log("Computing statistics (Pearson, Spearman, Linear Regression, Cosine, MSE/MAE) ...")
    pearson_corr, pearson_p = pearsonr(values1, values2)
    spearman_corr, spearman_p = spearmanr(values1, values2)
    slope, intercept, r_value, p_value, std_err = linregress(values1, values2)
    r_squared = r_value ** 2
    mse = mean_squared_error(values1, values2)
    mae = mean_absolute_error(values1, values2)
    cosine_sim = cosine_similarity(values1.values.reshape(1, -1),
                                   values2.values.reshape(1, -1))[0][0]

    log("==== Metrics ====")
    log(f"N         = {n_pairs}")
    log(f"R^2       = {r_squared:.6f}")
    log(f"Pearson   = {pearson_corr:.6f} (p={pearson_p:.3e})")
    log(f"Spearman  = {spearman_corr:.6f} (p={spearman_p:.3e})")
    log(f"Cosine    = {cosine_sim:.6f}")
    log(f"MSE       = {mse:.6e}")
    log(f"MAE       = {mae:.6e}")
    log(f"Trend     = y = {slope:.6f} * x + {intercept:.6f}")
    log("=================")

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(4, 4))

    # 散点（栅格化）
    sc = ax.scatter(
        values1, values2,
        alpha=args.alpha,
        marker=args.marker,
        s=16,
        linewidths=0.8,
        rasterized=True,
        color='#3dbfaa'
    )

    # 拟合线数据
    x_sorted = np.sort(values1.values)
    y_pred = slope * x_sorted + intercept

    # CI 阴影
    if args.ci:
        ci_level = args.ci / 100.0
        alpha_level = 1 - ci_level
        n = len(values1)
        y_err = values2 - (slope * values1 + intercept)
        residual_std_error = np.sqrt(np.sum(y_err**2) / (n - 2))
        t_val = t.ppf(1 - alpha_level / 2, df=n - 2)
        se_fit = residual_std_error * np.sqrt(1/n + (x_sorted - np.mean(values1))**2 /
                                              np.sum((values1 - np.mean(values1))**2))
        ci_upper = y_pred + t_val * se_fit
        ci_lower = y_pred - t_val * se_fit
        ax.fill_between(
            x_sorted, ci_lower, ci_upper,
            color='#1d8f8d', alpha=0.2, zorder=1,
            rasterized=True, label=f"{int(args.ci)}% CI"
        )

    # 趋势线（栅格化）
    ax.plot(
        x_sorted, y_pred,
        color='#1d8f8d', linewidth=1.5, label="Trend line",
        zorder=2
    )

    # x=y 参考线（矢量）
    vmin = float(min(values1.min(), values2.min()))
    vmax = float(max(values1.max(), values2.max()))
    ax.plot([vmin, vmax], [vmin, vmax],
            color="black", linestyle="--", linewidth=1.0, label="x = y")

    # 主图例
    main_legend = ax.legend(loc="upper left", frameon=True, title=None, fontsize=8)
    ax.add_artist(main_legend)

    # Metrics 文本
    metrics_lines = [
        f"{'R²':<9} = {r_squared:.3f}",
        f"{'Pearson':<9} = {pearson_corr:.3f}",
        f"{'Spearman':<9} = {spearman_corr:.3f}",
        f"{'Cosine':<9} = {cosine_sim:.3f}",
#       f"MSE       = {mse:.3e}",
#       f"MAE       = {mae:.3e}"
    ]
    metrics_text = "\n".join(metrics_lines)
    dummy_handle = Line2D([0], [0], color="white", alpha=0.0)

    # 自动位置选择
    loc = args.stats_loc if args.stats_loc else auto_metrics_location(values1, values2)

    if loc == "upper left":
        # 避开主图例（下移）
        second_legend = ax.legend(
            [dummy_handle], [metrics_text],
            loc="upper left", bbox_to_anchor=(0, 0.85),
            frameon=True, borderpad=1.0, handlelength=0,
            handletextpad=0.0, labelspacing=1.2,
            fontsize=8
        )
    else:
        second_legend = ax.legend(
            [dummy_handle], [metrics_text],
            loc=loc, frameon=True,
            borderpad=1.0, handlelength=0,
            handletextpad=0.0, labelspacing=1.2,
            fontsize=8
        )

    second_legend.get_frame().set_alpha(0.5)
    ax.add_artist(second_legend)

#   ax.set_title("Distance Concordance")
    ax.set_xlabel(args.xlabel)
    ax.set_ylabel(args.ylabel)
    ax.grid(True, linewidth=0.3, alpha=0.4)

    # 保存
    out_path = infer_output_path(args.output, args.format.lower())
    save_kwargs = {}
    if Path(out_path).suffix.lower() in {".png", ".tif", ".tiff", ".jpg", ".jpeg"}:
        save_kwargs["dpi"] = args.dpi
    plt.tight_layout()
    log(f"Saving figure to: {out_path}")
    plt.savefig(out_path, format=args.format.lower(), **save_kwargs)
    log(f"Saved figure to: {out_path}")

if __name__ == "__main__":
    main()

