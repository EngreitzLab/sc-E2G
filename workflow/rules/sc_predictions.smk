
# merge features with crispr data
rule overlap_features_crispr_apply:
	input:
		prediction_file = os.path.join(RESULTS_DIR, "{cluster}", "{model_name}", "encode_e2g_predictions.tsv.gz"),
		crispr = config['crispr_dataset'],
		feature_table_file = os.path.join(RESULTS_DIR, "{cluster}", "feature_table.tsv"),
		tss = encode_e2g_config['gene_TSS500']
	output: 
		features = os.path.join(RESULTS_DIR, "{cluster}", "{model_name}", "EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_{nafill}.tsv.gz"),
	conda:
		"../envs/sc_e2g.yml" 
	resources:
		mem_mb=32*1000
	script:
		"../scripts/merge_features_with_crispr_data_apply.R"

# calculate performance metrics 
rule crispr_benchmarking:
	input:
		crispr_features = expand(os.path.join(RESULTS_DIR, "{cluster}", "{model_name}", "EPCrisprBenchmark_ensemble_data_GRCh38.K562_features_NAfilled.tsv.gz"), zip, cluster=BIOSAMPLE_DF["cluster"], model_name=BIOSAMPLE_DF["model_dir_base"])
	output:
		comp_table = os.path.join(RESULTS_DIR, "crispr_benchmarking_performance_summary.tsv")
	params:
		model_names = BIOSAMPLE_DF["model_dir_base"].tolist(),
		model_thresh = BIOSAMPLE_DF["model_threshold"].tolist(),
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/sc_e2g.yml" 
	resources:
		mem_mb=64*1000
	shell: 
		""" 
		python {params.scripts_dir}/benchmark_performance.py \
			--crispr_features "{input.crispr_features}" \
			--output_file {output.comp_table} \
			--model_thresholds "{params.model_thresh}" \
			--model_names "{params.model_names}"
		"""

rule run_e2g_qnorm:
	input:
		final_features = os.path.join(RESULTS_DIR, "{cluster}", "genomewide_features.tsv.gz"),
	params:
		epsilon = 0.01,
		feature_table_file = lambda wildcards: encode_e2g.get_feature_table_file(wildcards.cluster, wildcards.model_name),
		trained_model = lambda wildcards: encode_e2g.get_trained_model(wildcards.cluster, wildcards.model_name),
		model_dir = lambda wildcards: encode_e2g._get_model_dir_from_wildcards(wildcards.cluster, wildcards.model_name, BIOSAMPLE_DF),
		crispr_benchmarking = config["benchmark_performance"],
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/sc_e2g.yml"
	resources:
		mem_mb=determine_mem_mb
	output: 
		prediction_file = os.path.join(RESULTS_DIR, "{cluster}", "{model_name}", "encode_e2g_predictions.tsv.gz")
	shell: 
		""" 
		python {params.scripts_dir}/run_e2g_cv.py \
			--predictions {input.final_features} \
			--feature_table_file {params.feature_table_file} \
			--epsilon {params.epsilon} \
			--trained_model {params.trained_model} \
			--model_dir {params.model_dir} \
			--crispr_benchmarking {params.crispr_benchmarking} \
			--output_file {output.prediction_file}
		"""