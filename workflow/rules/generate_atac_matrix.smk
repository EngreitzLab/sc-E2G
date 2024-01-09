## generate single-cell atac-seq matrix
rule generate_atac_matrix:
	input:
		kendall_pairs_path = 
			os.path.join(
				config["results_dir"], 
				"{cluster}", 
				"Kendall", 
				"Paires.tsv.gz"
			),
		multiome_sample_path = config["multiome_samples"],
		meta_data_path = config["meta_data"]
	output:
		atac_matrix = 
			os.path.join(
				config["results_dir"], 
				"{cluster}", 
				"Kendall", 
				"atac_matrix.csv.gz"
			)
	params:
		meta_col_cell_name = config["meta_data_col"]["cell_name"],
		meta_col_sample = config["meta_data_col"]["sample"],
		meta_col_cluster = config["meta_data_col"]["cluster"],
		meta_col_barcode = config["meta_data_col"]["barcode"]

	conda:
		"../envs/sc_e2g.yml"
	script:
		"../scripts/generate_atac_matrix.R"
