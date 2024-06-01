import pickle
import os
import click
import numpy as np
import pandas as pd
import scipy
from training_functions import statistic_aupr, statistic_precision_at_threshold, statistic_recall_at_threshold, threshold_at_target_recall


def performance_summary(cluster, model_name, model_threshold, crispr_features, n_boot=1000):
	# pct_missing
	if model_name.startswith("multiome"):
		n_zero = crispr_features.loc[crispr_features['ARC.E2G.Score']==0].shape[0]
	elif model_name.startswith("scATAC"):
		n_zero = crispr_features.loc[crispr_features['ABC.Score']==0].shape[0]
	pct_missing = n_zero/crispr_features.shape[0]
	
	# get scores
	Y_true = crispr_features['Regulated'].values.astype(np.int64)
	#Y_pred = crispr_features['ENCODE-rE2G.Score.cv.qnorm']
	Y_pred = crispr_features['ENCODE-rE2G.Score.qnorm']


	# auprc
	res_aupr =  scipy.stats.bootstrap((Y_true, Y_pred), statistic_aupr, n_resamples=n_boot, paired=True, confidence_level=0.95, method='BCa')

	# precision at 70% recall
	thresh_70 = threshold_at_target_recall(Y_true, Y_pred, 0.7) # will return None if max recall < 70%
	if thresh_70 is not None:
		res_70_prec = scipy.stats.bootstrap((Y_true, Y_pred),  lambda Y_true, Y_pred: statistic_precision_at_threshold(Y_true, Y_pred, thresh_70), n_resamples=n_boot, paired=True, confidence_level=0.95, method='BCa')
	prec_70_mean = 0 if thresh_70 is None else np.mean(res_70_prec.bootstrap_distribution)
	prec_70_low = 0 if thresh_70 is None else res_70_prec.confidence_interval[0]
	prec_70_high = 0 if thresh_70 is None else res_70_prec.confidence_interval[1]

	# precision at 50% recall
	thresh_50 = threshold_at_target_recall(Y_true, Y_pred, 0.5) # will return None if max recall < 70%
	if thresh_50 is not None:
		res_50_prec = scipy.stats.bootstrap((Y_true, Y_pred),  lambda Y_true, Y_pred: statistic_precision_at_threshold(Y_true, Y_pred, thresh_50), n_resamples=n_boot, paired=True, confidence_level=0.95, method='BCa')
	prec_50_mean = 0 if thresh_50 is None else np.mean(res_50_prec.bootstrap_distribution)
	prec_50_low = 0 if thresh_50 is None else res_50_prec.confidence_interval[0]
	prec_50_high = 0 if thresh_50 is None else res_50_prec.confidence_interval[1]

	# precision at model threshold
	model_thresh = float(model_threshold)
	res_prec_model_thresh = scipy.stats.bootstrap((Y_true, Y_pred),  lambda Y_true, Y_pred: statistic_precision_at_threshold(Y_true, Y_pred, model_thresh), n_resamples=n_boot, paired=True, confidence_level=0.95, method='BCa')
	prec_mean_model_thresh =  np.mean(res_prec_model_thresh.bootstrap_distribution)
	prec_low_model_thresh =  res_prec_model_thresh.confidence_interval[0]
	prec_high_model_thresh =  res_prec_model_thresh.confidence_interval[1]

	# recall at model threshold
	res_recall_model_thresh = scipy.stats.bootstrap((Y_true, Y_pred),  lambda Y_true, Y_pred: statistic_recall_at_threshold(Y_true, Y_pred, model_thresh), n_resamples=n_boot, paired=True, confidence_level=0.95, method='BCa')
	recall_mean_model_thresh =  np.mean(res_recall_model_thresh.bootstrap_distribution)
	recall_low_model_thresh =  res_recall_model_thresh.confidence_interval[0]
	recall_high_model_thresh =  res_recall_model_thresh.confidence_interval[1]

	res_row = pd.DataFrame({'cluster': cluster, 'model': model_name,  'AUPRC': np.mean(res_aupr.bootstrap_distribution), 'AUPRC_95CI_low': res_aupr.confidence_interval[0], 'AUPRC_95CI_high': res_aupr.confidence_interval[1],
								 'precision_70_pct_recall': prec_70_mean, 'precision_70_pct_recall_95CI_low': prec_70_low, 'precision_70_pct_recall_95CI_high': prec_70_high, 'threshold_70_pct_recall': thresh_70,
								 'precision_50_pct_recall': prec_50_mean, 'precision_50_pct_recall_95CI_low': prec_50_low, 'precision_50_pct_recall_95CI_high': prec_50_high, 'threshold_50_pct_recall': thresh_50,
								 'precision_model_threshold': prec_mean_model_thresh, 'precision_model_threshold_95CI_low': prec_low_model_thresh, 'precision_model_threshold_95CI_high': prec_high_model_thresh,
								 'recall_model_threshold': recall_mean_model_thresh, 'recall_model_threshold_95CI_low': recall_low_model_thresh, 'recall_model_threshold_95CI_high': recall_high_model_thresh,
								 'pct_missing_elements': pct_missing}, index=[0])
	return res_row

	
@click.command()
@click.option("--crispr_features", required=True)
@click.option("--output_file", required=True)
@click.option("--model_names", required=True)
@click.option("--model_thresholds", required=True)

def main(crispr_features, output_file, model_names, model_thresholds):
	# organize inputs
	crispr_features = crispr_features.split(" ") # .../cluster/model_name/crispr_features.tsv.gz
	clusters = [file_path.split(os.sep)[-3] for file_path in crispr_features]
	models = [file_path.split(os.sep)[-2] for file_path in crispr_features]
	preds = pd.DataFrame({'cluster': clusters, 'model_name': models, 'crispr_file': crispr_features})

	model_dict = {model_name: float(threshold) for model_name, threshold in zip(model_names.split(" "), model_thresholds.split(" "))}
	preds['model_threshold'] = [model_dict[name] for name in preds['model_name']]

	# initiate final df
	df = pd.DataFrame(columns = ['cluster', 'model', 'AUPRC', 'AUPRC_95CI_low', 'AUPRC_95CI_high',
								 'precision_70_pct_recall', 'precision_70_pct_recall_95CI_low', 'precision_70_pct_recall_95CI_high', 'threshold_70_pct_recall', 
								 'precision_50_pct_recall', 'precision_50_pct_recall_95CI_low', 'precision_50_pct_recall_95CI_high', 'threshold_50_pct_recall', 
								 'precision_model_threshold', 'precision_model_threshold_95CI_low', 'precision_model_threshold_95CI_high',
								 'recall_model_threshold', 'recall_model_threshold_95CI_low', 'recall_model_threshold_95CI_high',
								 'pct_missing_elements'])
	
	# iterate through list
	for row in preds.itertuples(index=False):
		print(row)
		crispr_features = pd.read_csv(row.crispr_file, sep="\t")
		res_row = performance_summary(row.cluster, row.model_name, row.model_threshold, crispr_features)
		df = pd.concat([df, res_row])
			
	df = df.sort_values(by='AUPRC', ascending=False)
	df.to_csv(output_file, sep = '\t', index=False)
	
if __name__ == "__main__":
	main()
