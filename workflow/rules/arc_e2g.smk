## integrate Kendall and ABC into ARC-E2G
rule arc_e2g:
	input:
		abc_predictions = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"Predictions", 
				"EnhancerPredictionsAllPutative.tsv.gz"
			),
		kendall_predictions = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"Kendall", 
				"Paires.Kendall.tsv.gz")
	output:
		arc_predictions = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"ARC",
				"EnhancerPredictionsAllPutative_ARC.tsv.gz")
	conda:
		"../envs/sc_e2g.yml"
	script:
		"../scripts/compute_arc_e2g.R"
