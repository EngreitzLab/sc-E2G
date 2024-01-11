## Convert fragment file to tagAlign file
rule frag_to_tagAlign:
	input:
		frag_file = 
			os.path.join(
				config["results_dir"], 
				"{cluster}", 
				"Fragments", 
				"atac_fragments.tsv.gz"
			)
	params:
		chrom_sizes = encode_e2g.ABC.config['ref']['chrom_sizes']
	output:
		tagAlign_sort_file = 
			os.path.join(
				config["results_dir"], 
				"{cluster}", 
				"Fragments",
				"tagAlign.sort.gz"
			)
	conda:
		"../envs/sc_e2g.yml"
	shell:
		"""
		# Make, sort and compress tagAlign file from fragment file
		gunzip -c {input.frag_file} | \
		awk 'BEGIN{{OFS="\\t"}} {{
			mid = int(($2 + $3) / 2); # Calculate middle point
			# Repeat the output tags according to the 5th column of the fragment file
			for(i = 1; i <= $5; i++) {{
				print $1, $2, mid, "N", 1000, "+"; # positive strand region
				print $1, mid + 1, $3, "N", 1000, "-"; # negative strand region
			}}
		}}' | \
		bedtools sort -faidx {params.chrom_sizes} -i stdin | \
		bgzip > {output.tagAlign_sort_file}
		
		# Index the tagAlign file
		tabix -p bed {output.tagAlign_sort_file}
		""" 

