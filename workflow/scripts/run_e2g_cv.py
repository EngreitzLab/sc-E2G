import pickle
import os
import click
import numpy as np
import pandas as pd
from scipy import interpolate

MODEL = "ENCODE-rE2G"

def make_e2g_predictions(df_enhancers, feature_list, trained_model, epsilon):
	# transform the features
	X = df_enhancers.loc[:, feature_list]
	X = np.log(np.abs(X) + epsilon)

	with open(trained_model, "rb") as f:
		model = pickle.load(f)
	probs = model.predict_proba(X)
	df_enhancers[MODEL + ".Score"] = probs[:, 1]
	return df_enhancers
	
def make_e2g_predictions_cv(df_enhancers, feature_list, cv_models, epsilon):
	score_name = MODEL + '.Score.cv'
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
			df_enhancers.loc[idx_test, score_name] = probs[:,1]

	return df_enhancers

def qnorm_scores(df_enhancers, qnorm_ref, crispr_benchmarking):
	# generate reference
	ref_quantiles = np.linspace(0.01, 0.99, 99, endpoint=True) # [0.01, 0.02, ... , 0.99]
	ref_scores = np.quantile(qnorm_ref, ref_quantiles)
	interpfunc = interpolate.interp1d(ref_quantiles, ref_scores, kind="linear", fill_value="extrapolate")

	# qnorm regular score
	score_quantile = df_enhancers[MODEL + '.Score'].rank() / float(len(df_enhancers))
	df_enhancers[MODEL + '.Score.qnorm'] = interpfunc(score_quantile).clip(0)
	print(df_enhancers[MODEL + '.Score'].describe())
	print(df_enhancers[MODEL + '.Score.qnorm'].describe())

	# qnorm cv score
	if crispr_benchmarking:
		cv_score_quantile = df_enhancers[MODEL + '.Score.cv'].rank() / float(len(df_enhancers))
		df_enhancers[MODEL + '.Score.cv.qnorm'] = interpfunc(cv_score_quantile).clip(0)

	return df_enhancers
	
@click.command()
@click.option("--predictions", required=True)
@click.option("--feature_table_file", required=True)
@click.option("--trained_model", required=True)
@click.option("--model_dir", required=True)
@click.option("--crispr_benchmarking", required=True)
@click.option("--epsilon", type=float, default=0.01)
@click.option("--output_file", required=True)
def main(predictions, feature_table_file, trained_model, model_dir, crispr_benchmarking, epsilon, output_file):
	feature_table = pd.read_csv(feature_table_file, sep="\t")
	feature_list = feature_table["feature"]

	# genomewide features
	df_enhancers = pd.read_csv(predictions, sep="\t") 
	df_enhancers = df_enhancers.replace([np.inf, -np.inf], np.nan)
	df_enhancers = df_enhancers.fillna(0)

	# read in qnorm ref
	qnorm_file = os.path.join(model_dir, "qnorm_reference.tsv.gz") # col name = "ENCODE-rE2G.Score'"
	qnorm_ref = pd.read_csv(qnorm_file, sep="\t")
	qnorm_ref = qnorm_ref['ENCODE-rE2G.Score'].to_numpy()

	df_enhancers = make_e2g_predictions(df_enhancers, feature_list, trained_model, epsilon)  # add "ENCODE-rE2G.Score"
	
	if crispr_benchmarking:
		cv_models = os.path.join(model_dir, "cv_models")
		df_enhancers = make_e2g_predictions_cv(df_enhancers, feature_list, cv_models, epsilon) # add "ENCODE-rE2G.Score.cv"

	df_enhancers = qnorm_scores(df_enhancers, qnorm_ref, crispr_benchmarking) # add "ENCODE-rE2G.Score.qnorm" &  if crispr_benchmarking, "ENCODE-rE2G.Score.cv.qnorm"

	df_enhancers.to_csv(output_file, compression="gzip", sep="\t", index=False)

if __name__ == "__main__":
	main()
