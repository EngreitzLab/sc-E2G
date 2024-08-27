import pickle
import os
import click
import numpy as np
import pandas as pd
from scipy import interpolate

SCORE_COL_BASE = "E2G.Score"

def make_e2g_predictions(df_enhancers, feature_list, trained_model, tpm_threshold, epsilon):
	# transform the features
	X = df_enhancers.loc[:, feature_list]
	print(df_enhancers.head())
	print(X.head())
	X = np.log(np.abs(X) + epsilon)

	with open(trained_model, "rb") as f:
		model = pickle.load(f)
	probs = model.predict_proba(X)

	if ("RNA_pseudobulkTPM" in df_enhancers.columns) and (tpm_threshold > 0):
		df_enhancers[SCORE_COL_BASE  + ".ignoreTPM"] = probs[:, 1]
		df_enhancers[SCORE_COL_BASE] = [score if tpm >= tpm_threshold else 0 
			for (score, tpm) in zip(df_enhancers[SCORE_COL_BASE + ".ignoreTPM"] , df_enhancers["RNA_pseudobulkTPM"])]
	else:
		df_enhancers[SCORE_COL_BASE] = probs[:, 1]

	return df_enhancers
	
def make_e2g_predictions_cv(df_enhancers, feature_list, cv_models, tpm_threshold, epsilon):
	score_col = SCORE_COL_BASE + '.cv'
	tpm_filter = ("RNA_pseudobulkTPM" in df_enhancers.columns) and (tpm_threshold > 0)

	# get scores per chrom
	X = df_enhancers.loc[:,feature_list]
	X = np.log(np.abs(X) + epsilon)
	chr_list = np.unique(df_enhancers['chr'])

	for chr in chr_list:
		idx_test = df_enhancers[df_enhancers['chr']==chr].index.values
		if len(idx_test)>0:
			X_test = X.loc[idx_test,:]
			with open(os.path.join(cv_models, f"model_test_{chr}.pkl"), 'rb') as f:
				model = pickle.load(f)
			probs = model.predict_proba(X_test)
			if tpm_filter:
				df_enhancers[idx_test, score_col  + ".ignoreTPM"] = probs[:, 1]
				df_enhancers[idx_test, score_col] = [score if tpm >= tpm_threshold else 0 
					for (score, tpm) in zip(df_enhancers.loc[idx_test, score_col + ".ignoreTPM"] , df_enhancers.loc[idx_test, "RNA_pseudobulkTPM"])]
			else:
				df_enhancers.loc[idx_test, score_col] = probs[:,1]

	return df_enhancers


def calculate_quantiles(score_column):
	score_column_zero_replaced = score_column.replace(0, np.nan)
	ranks = score_column_zero_replaced.rank(method='average', na_option='top') # Calculate the rank, treating NaNs as the lowest possible ranks
	quantiles = (ranks - 1) / (len(score_column) - 1) # Calculate the quantiles
	quantiles[score_column == 0] = 0 # Replace NaN ranks (original 0 values) with 0th quantile
	return quantiles

def qnorm_scores(df_enhancers, qnorm_ref, crispr_benchmarking):
	# generate reference
	#ref_quantiles = np.linspace(0.01, 0.99, 99, endpoint=True) # [0.01, 0.02, ... , 0.99]
	#n_unique_scores = len(np.unique(qnorm_ref))
	n_unique_scores = 100000
	print(n_unique_scores)
	ref_quantiles = np.linspace(0, 1, n_unique_scores, endpoint=True) # [0, 0.01, 0.02, ... , 0.99, 1]
	ref_scores = np.quantile(qnorm_ref, ref_quantiles)
	interpfunc = interpolate.interp1d(ref_quantiles, ref_scores, kind="linear", fill_value="extrapolate")

	# qnorm regular score
	#score_quantile = df_enhancers[SCORE_COL_BASE + '.Score'].rank() / float(len(df_enhancers))
	score_quantile = calculate_quantiles(df_enhancers[SCORE_COL_BASE])
	df_enhancers[SCORE_COL_BASE + '.qnorm'] = interpfunc(score_quantile).clip(0)
	print(df_enhancers[SCORE_COL_BASE].describe())
	print(df_enhancers[SCORE_COL_BASE + '.qnorm'].describe())

	# qnorm cv score
	if crispr_benchmarking:
		cv_score_quantile = calculate_quantiles(df_enhancers[SCORE_COL_BASE + '.cv'])
		df_enhancers[SCORE_COL_BASE + '.cv.qnorm'] = interpfunc(cv_score_quantile).clip(0)

	return df_enhancers
	
@click.command()
@click.option("--predictions", required=True)
@click.option("--feature_table_file", required=True)
@click.option("--trained_model", required=True)
@click.option("--crispr_benchmarking", required=True)
@click.option("--epsilon", type=float, default=0.01)
@click.option("--tpm_threshold", type=float, default=0)
@click.option("--qnorm_file", required=True)


def main(predictions, feature_table_file, trained_model, crispr_benchmarking, epsilon, tpm_threshold, qnorm_file):
	feature_table = pd.read_csv(feature_table_file, sep="\t")
	feature_list = feature_table["feature"]

	# genomewide features
	df_enhancers = pd.read_csv(predictions, sep="\t") 
	df_enhancers = df_enhancers.replace([np.inf, -np.inf], np.nan)
	df_enhancers = df_enhancers.fillna(0)
	df_enhancers = df_enhancers.replace({True: 1, False: 0})

	df_enhancers = make_e2g_predictions(df_enhancers, feature_list, trained_model, tpm_threshold, epsilon)  # add "E2G.Score" and potentially "E2G.Score.ignoreTPM"
	
	# if crispr_benchmarking:
	# 	cv_models = os.path.join(model_dir, "cv_models")
	# 	df_enhancers = make_e2g_predictions_cv(df_enhancers, feature_list, cv_models, tpm_threshold, epsilon) # add "E2G.Score.cv" and potentially "E2G.Score.cv.gnoreTPM"

	#df_enhancers = qnorm_scores(df_enhancers, qnorm_ref, crispr_benchmarking) # add "E2G.Score.qnorm" &  if crispr_benchmarking, "E2G.Score.cv.qnorm"

	df_enhancers["E2G.Score"].to_csv(qnorm_file, compression="gzip", sep="\t", index=False)

if __name__ == "__main__":
	main()
