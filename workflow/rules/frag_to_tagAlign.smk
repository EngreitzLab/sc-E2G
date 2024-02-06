## Create link from fragment file to tagAlign file
rule frag_to_tagAlign:
	input:
		frag_file = 
			lambda wildcards: cell_cluster_config.loc[wildcards.cluster, "atac_frag_file"],
		frag_file_index = 
			lambda wildcards: cell_cluster_config.loc[wildcards.cluster, "atac_frag_file"] + ".tbi"
	output:
		tagAlign_sort_file = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"tagAlign",
				"tagAlign.sort.gz"
			),
		tagAlign_sort_file_index = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"tagAlign", 
				"tagAlign.sort.gz.tbi"
			)
	shell:
		"""
		ln {input.frag_file} {output.tagAlign_sort_file}
		ln {input.frag_file_index} {output.tagAlign_sort_file_index}
		""" 

