## Add external features for training
rule add_external_features_training:
	input:
		predictions_extended = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}",
				"{model}", 
				"ActivityOnly_features.tsv.gz"
			),
		arc_predictions = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}",
				"ARC", 
				"EnhancerPredictionsAllPutative_ARC.tsv.gz"
			),
		feature_table_file = 
			lambda wildcards: model_config.loc[wildcards.model, "feature_table"]
	output:
		plus_external_features = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}",
				"{model}", 
				"ActivityOnly_plus_external_features.tsv.gz")
	conda:
		"../envs/sc_e2g.yml"
	script:
		"../scripts/integrate_kendall.R"
