## compute Kendall correlation
rule compute_kendall:
	input:
		kendall_pairs_path = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"Kendall", 
				"Pairs.tsv.gz"
			),
		atac_matix = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"Kendall", 
				"atac_matrix.csv.gz"
			),
		rna_matix = 
			lambda wildcards: CELL_CLUSTER_DF.loc[wildcards.cluster, "rna_matrix_file"]
	output:
		kendall_predictions = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"Kendall", 
				"Pairs.Kendall.tsv.gz") 
	resources: 
		mem_mb=200*1000
	conda:
		"../envs/sc_e2g.yml"
	script:
		"../scripts/compute_kendall.R"
