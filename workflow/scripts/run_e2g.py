import pickle

import click
import numpy as np
import pandas as pd

MODEL = "ENCODE-E2G"


def make_e2g_predictions(df_enhancers, feature_list, models_dir, epsilon):
    # transform the features
    X = df_enhancers.loc[:, feature_list]
    X.columns = feature_list
    X = np.log(np.abs(X) + epsilon)
    chr_list = np.unique(df_enhancers["chr"])

    for chr in chr_list:
        idx_test = df_enhancers[df_enhancers["chr"] == chr].index.values

        if len(idx_test) > 0:
            X_test = X.loc[idx_test, :]

            with open(f"{models_dir}/model_{MODEL}_test_{chr}.pkl", "rb") as f:
                model = pickle.load(f)

            probs = model.predict_proba(X_test)
            df_enhancers.loc[idx_test, MODEL + ".Score"] = probs[:, 1]
    return df_enhancers


@click.command()
@click.option("--predictions", required=True)
@click.option("--feature_table_file", required=True)
@click.option("--models_dir", required=True)
@click.option("--epsilon", type=float, default=0.01)
@click.option("--output_file", required=True)
def main(predictions, feature_table_file, models_dir, epsilon, output_file):
    feature_table = pd.read_csv(feature_table_file, sep="\t")
    feature_list = feature_table[feature_table[MODEL] == 1]["features"]

    df_enhancers = pd.read_csv(predictions, sep="\t")
    df_enhancers = df_enhancers.replace([np.inf, -np.inf], np.nan)
    df_enhancers = df_enhancers.fillna(0)

    df_enhancers = make_e2g_predictions(df_enhancers, feature_list, models_dir, epsilon)
    df_enhancers.to_csv(output_file, compression="gzip", sep="\t", index=False)


if __name__ == "__main__":
    main()
