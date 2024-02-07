## Convert fragment file to tagAlign file
rule frag_to_tagAlign:
	input:
		frag_file = 
			lambda wildcards: cell_cluster_config.loc[wildcards.cluster, "atac_frag_file"]
	output:
		tagAlign_sort_file = 
			os.path.join(
				RESULTS_DIR, 
				"{cluster}", 
				"tagAlign",
				"tagAlign.sort.gz"
			)
	conda:
		"../envs/sc_e2g.yml"
	threads: 5
	shell:
		"""
		# Make, sort and compress tagAlign file from fragment file
		LC_ALL=C zcat {input.frag_file}  | sed '/^#/d' | \
		awk -v OFS='\t' '{{mid=int(($2+$3)/2); print $1,$2,mid,"N",1000,"+"; print $1,mid+1,$3,"N",1000,"-"}}' | \
		sort -k 1,1V -k 2,2n -k3,3n --parallel 5 | \
		bgzip -c > {output.tagAlign_sort_file}  

		# Index the tagAlign file
		tabix -p bed {output.tagAlign_sort_file}
		""" 
