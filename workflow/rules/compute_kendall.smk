from functools import partial

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
				"atac_matrix.rds"
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
		mem_mb=partial(determine_mem_mb, min_gb=63),
		#mem_mb=126000,
		runtime=12*60
	conda:
		"../envs/sc_e2g.yml"
	script:
		"../scripts/compute_kendall.R"
