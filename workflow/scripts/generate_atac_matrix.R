## Generate single-cell ATAC-seq matrix for a specific cluster

# Load required packages
suppressPackageStartupMessages({
  library(genomation)
  library(GenomicRanges)
  library(Signac)
})

# Import parameters from Snakemake
kendall_pairs_path = snakemake@input[["kendall_pairs_path"]]
multiome_sample_path = snakemake@input[["multiome_sample_path"]]
meta_data_path = snakemake@input[["meta_data_path"]]
atac_matrix_path = snakemake@output[["atac_matrix"]]
meta_col_cell_name = snakemake@params[["meta_col_cell_name"]]
meta_col_sample = snakemake@params[["meta_col_sample"]]
meta_col_cluster = snakemake@params[["meta_col_cluster"]]
meta_col_barcode = snakemake@params[["meta_col_barcode"]]
cluster = snakemake@wildcards$cluster

# Read the enhancer-gene pairs and extract unique peaks
pairs.e2g = readGeneric(kendall_pairs_path,
                        keep.all.metadata = T,
                        header = T)
bed.peaks = pairs.e2g[!duplicated(mcols(pairs.e2g)[,"PeakName"])]
mcols(bed.peaks) = NULL

# Read the multiome sample information and meta_data
df.multiome_sample = read.delim(multiome_sample_path)
df.meta_data = read.csv(meta_data_path)

# Create a list to store Signac Fragment objects for each sample
list.fragments = list()
for (n in 1:nrow(df.multiome_sample)) {
  # Extract sample-specific metadata
  sample.tmp = df.multiome_sample[n,"sample"]
  df.meta_data.tmp = 
    df.meta_data[df.meta_data[,meta_col_sample] == sample.tmp,]
  
  # Assign cell names to barcodes
  cells.tmp = df.meta_data.tmp[,meta_col_barcode]
  names(cells.tmp) = df.meta_data.tmp[,meta_col_cell_name]
  
  # Create a Signac Fragment object for each sample
  list.fragments[[n]] =
    CreateFragmentObject(path = df.multiome_sample[n,"atac_frag_path"],
                         cells = cells.tmp)
}

# Construct the ATAC-seq matrix for the specified cluster
atac.matrix <- FeatureMatrix(
  fragments = list.fragments,
  features = bed.peaks,
  cells = 
    df.meta_data[df.meta_data[,meta_col_cluster] == cluster,
                 meta_col_cell_name]
)

# Write the ATAC-seq matrix to a gzipped file
write.csv(as.matrix(atac.matrix),
          file = gzfile(atac_matrix_path),
          quote = F)

