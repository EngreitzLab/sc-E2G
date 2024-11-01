import pickle
import os
import click
import numpy as np
import pandas as pd
from scipy import interpolate


def get_qnorm_quantiles(input_df, n_scores):
    all_scores = input_df["E2G.Score"]
    n_unique_scores = int(n_scores)

    ref_quantiles = np.linspace(
        0, 1, n_unique_scores, endpoint=True
    )  # [0, 0.01, 0.02, ... , 0.99, 1]
    ref_scores = np.quantile(all_scores, ref_quantiles)

    out_df = pd.DataFrame({"quantile": ref_quantiles, "reference_score": ref_scores})

    return out_df


@click.command()
@click.option("--input_file", required=True)
@click.option("--output_file", required=True)
@click.option("--n_scores", type=float, required=True)
def main(input_file, output_file, n_scores):
    input_df = pd.read_csv(input_file, sep="\t")
    out_df = get_qnorm_quantiles(input_df, n_scores)
    out_df.to_csv(output_file, compression="gzip", sep="\t", index=False)


if __name__ == "__main__":
    main()
