## generate single-cell atac-seq matrix
rule generate_atac_matrix:
	input:
		kendall_pairs_path = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"Kendall", 
				"Paires.tsv.gz"
			),
		atac_frag_path = 
			lambda wildcards: cell_cluster_config.loc[wildcards.cluster, "atac_frag_file"],
		rna_matrix_path = 
			lambda wildcards: cell_cluster_config.loc[wildcards.cluster, "rna_matrix_file"]
	output:
		atac_matrix_path = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"Kendall", 
				"atac_matrix.csv.gz"
			)
	conda:
		"../envs/sc_e2g.yml"
	script:
		"../scripts/generate_atac_matrix.R"
