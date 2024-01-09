## compute Kendall correlation
rule compute_kendall:
	input:
		kendall_pairs_path = 
			os.path.join(
				config["results_dir"], 
				"{cluster}", 
				"Kendall", 
				"Paires.tsv.gz"
			),
		atac_matix = 
			os.path.join(
				config["results_dir"], 
				"{cluster}", 
				"Kendall", 
				"atac_matrix.csv.gz"
			),
		rna_matix = config["rna_matrix"]
	output:
		kendall_predictions = 
			os.path.join(
				config["results_dir"], 
				"{cluster}", 
				"Kendall", 
				"Paires.Kendall.tsv.gz") 
	conda:
		"../envs/sc_e2g.yml"
	script:
		"../scripts/compute_kendall.R"
