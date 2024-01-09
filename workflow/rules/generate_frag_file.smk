## Generate fragment files for each cell cluster
rule generate_frag_file:
	input:
		multiome_sample_path = config["multiome_samples"],
		meta_data_path = config["meta_data"]
	output:
		frag_file = 
			os.path.join(
				config["results_dir"], 
				"{cluster}", 
				"Fragments", 
				"atac_fragments.tsv.gz"
			)
	params:
		meta_col_sample = config["meta_data_col"]["sample"],
		meta_col_cluster = config["meta_data_col"]["cluster"],
		meta_col_barcode = config["meta_data_col"]["barcode"]
	conda:
		"../envs/sc_e2g.yml"
	script:
		"../scripts/generate_frag_file.R"
