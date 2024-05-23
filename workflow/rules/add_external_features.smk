## Add external features
rule add_external_features:
	input:
		predictions_extended = 
			os.path.join(
				RESULTS_DIR, 
				"{biosample}", 
				"ActivityOnly_features.tsv.gz"
			),
		arc_predictions = 
			os.path.join(
				RESULTS_DIR, 
				"{biosample}", 
				"ARC", 
				"EnhancerPredictionsAllPutative_ARC.tsv.gz"
			),
		feature_table_file = 
			os.path.join(RESULTS_DIR, "{biosample}", "feature_table.tsv")
	output:
		plus_external_features = 
			os.path.join(
				RESULTS_DIR, 
				"{biosample}", 
				"ActivityOnly_plus_external_features.tsv.gz")
	conda:
		"../envs/sc_e2g.yml"
	resources:
		mem_mb=32*1000
	script:
		"../scripts/integrate_kendall.R"
