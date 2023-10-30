import click
import glob
import os
import pandas as pd
from typing import Dict
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

STATS_SUFFIX = "_stats.tsv"
VALUE_KEY = "Value"
METRICS = [
    "num_enh",
    "num_genes",
    "num_enh_gene_links",
    "num_genes_with_1_enh_min",
    "mean_num_genes_per_enh",
    "mean_num_enh_per_gene",
    "mean_num_enh_per_gene_no_prom",
    "mean_log10_dist_to_tss",
    "mean_enh_region_size"
]

def load_stat_files(stat_files) -> Dict[str, pd.DataFrame]:
    results = {}
    for stat_file in stat_files:
        cell_cluster = os.path.basename(os.path.dirname(stat_file))
        results[cell_cluster] = pd.read_csv(stat_file, sep="\t").set_index("Metric")
    return results

def snake_to_title(snake_case):
    words = snake_case.split('_')
    # Capitalize the first letter of each word except the first one
    title = ' '.join(word.capitalize() for word in words)
    return title

def plot_violin(stat_dfs, metric):
    plt.clf()
    points = [df.loc[metric, VALUE_KEY] for df in stat_dfs.values()]
    mean, median = np.mean(points), np.median(points)
    ax = sns.violinplot(y=points)
    ax.set_title(snake_to_title(metric))
    plt.scatter([], [], label=f"n={len(points)}\nMean={mean:.2f}\nMedian={median:.2f}")
    plt.legend()
    return ax.get_figure()

def plot_scatter(stat_dfs, metadata_df, metric):
    plt.clf()
    X, Y = [], []
    K562_X, K562_Y= [], []
    for cell_cluster in stat_dfs:
        row = metadata_df[metadata_df["CellClusterID"] == cell_cluster].iloc[0]
        num_fragments = row["nCells"] * row["MeanATACFragmentsPerCell"]
        metric_val = stat_dfs[cell_cluster].loc[metric, VALUE_KEY]
        if "K562" in row["ManualAnnotationLabel"]:
            K562_X.append(int(num_fragments))
            K562_Y.append(metric_val)
        else:
            X.append(int(num_fragments))
            Y.append(metric_val)
            

    ax = sns.scatterplot(x=X, y=Y, label=f"n={len(X)}")
    ax = sns.scatterplot(x=K562_X, y=K562_Y, color="red", label=f"K562 n={len(K562_X)}")
    plt.xscale("log")
    plt.axvline(x=2e6, color="red", linestyle="-", label=f"2 million fragments")
    ax.set_xlabel("Num Fragments")
    ax.set_ylabel(snake_to_title(metric))
    ax.set_title(snake_to_title(metric))
    plt.legend()

    return ax.get_figure()
    

@click.command()
@click.option("--results_dir", type=str, required=True)
@click.option("--output_file", type=str, default="qc_plots.pdf")
@click.option("--y2ave_metadata", type=str)
def main(results_dir, output_file, y2ave_metadata):
    if y2ave_metadata:
        metadata_df = pd.read_csv(y2ave_metadata, sep="\t")

    stat_files = glob.glob(os.path.join(results_dir, '*', f"*{STATS_SUFFIX}"))
    stat_dfs = load_stat_files(stat_files)
    with PdfPages(output_file) as pdf_writer:
        for metric in METRICS:
            pdf_writer.savefig(plot_violin(stat_dfs, metric))
            if y2ave_metadata:
                pdf_writer.savefig(plot_scatter(stat_dfs, metadata_df, metric))

if __name__ == "__main__":
    main()