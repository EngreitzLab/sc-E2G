## compute Kendall correlation
rule compute_kendall:
	input:
		kendall_pairs_path = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"Kendall", 
				"Paires.tsv.gz"
			),
		atac_matix = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"Kendall", 
				"atac_matrix.csv.gz"
			),
		rna_matix = 
			lambda wildcards: cell_cluster_config.loc[wildcards.cluster, "rna_matrix_file"]
	output:
		kendall_predictions = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"Kendall", 
				"Paires.Kendall.tsv.gz") 
	conda:
		"../envs/sc_e2g.yml"
	script:
		"../scripts/compute_kendall.R"
