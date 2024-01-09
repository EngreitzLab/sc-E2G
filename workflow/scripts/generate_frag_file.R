## generate fragment files for each cell cluster

# Import parameters from Snakemake
multiome_sample_path = snakemake@input[["multiome_sample_path"]]
meta_data_path = snakemake@input[["meta_data_path"]]
output_frag_file_path = snakemake@output[["frag_file"]]
meta_col_sample = snakemake@params[["meta_col_sample"]]
meta_col_cluster = snakemake@params[["meta_col_cluster"]]
meta_col_barcode = snakemake@params[["meta_col_barcode"]]
cluster = snakemake@wildcards$cluster

# Read the multiome sample information and meta_data
df.multiome_sample = read.delim(multiome_sample_path,
                                row.names = "sample")
df.meta_data = read.csv(meta_data_path)

# Initialize a list to store filtered fragments for each sample
list.frag = list()

for (sample in rownames(df.multiome_sample)) {
  # Filter the metadata based on the current sample and cluster
  df.meta_data.filter = 
    df.meta_data[df.meta_data[,meta_col_sample] == sample &
                   df.meta_data[,meta_col_cluster] == cluster, ]

  # Read fragments for the current sample
  df.frag = read.table(df.multiome_sample[sample,"atac_frag_path"])

  # Filter fragments based on barcodes matching the meta_data
  list.frag[[sample]] = 
    df.frag[df.frag[,4] %in% df.meta_data.filter[,meta_col_barcode],]
}

# Combine the filtered fragments from all samples
frag.output = do.call(rbind,list.frag)

# Sort the combined fragments based on coordinates
frag.output.sort = frag.output[order(frag.output[,2]),]
frag.output.sort = frag.output.sort[order(frag.output.sort[,1]),]

# Write the sorted fragments to a gzipped file
write.table(frag.output.sort,
            file = gzfile(output_frag_file_path),
            row.names = F,
            col.names = F,
            sep = "\t",
            quote = F)

