## Generate single-cell ATAC-seq matrix for a specific cluster

# Load required packages
suppressPackageStartupMessages({
  library(genomation)
  library(GenomicRanges)
  library(Signac)
  library(anndata)
  library(tools)
})

# Import parameters from Snakemake
kendall_pairs_path = snakemake@input[["kendall_pairs_path"]]
atac_frag_path = snakemake@input[["atac_frag_path"]]
rna_matrix_path = snakemake@input[["rna_matrix_path"]]
atac_matrix_path = snakemake@output[["atac_matrix_path"]]

# Read the enhancer-gene pairs and extract unique peaks
pairs.e2g = readGeneric(kendall_pairs_path,
                        keep.all.metadata = T,
                        header = T)
bed.peaks = pairs.e2g[!duplicated(mcols(pairs.e2g)[,"PeakName"])]
mcols(bed.peaks) = NULL

# Read the rna matrix to extract cell name
if (file_ext(rna_matrix_path) == "h5ad") {
  rna_matrix <- t(read_h5ad(rna_matrix_path)$X)
} else {
  rna_matrix = read.csv(rna_matrix_path,
                        row.names = 1,
                        check.names = F)
}

# Create a list to store Signac Fragment object
list.fragments = list()
cells.use = colnames(rna_matrix)
names(cells.use) = colnames(rna_matrix)
rm(rna_matrix)
list.fragments[[1]] =
  CreateFragmentObject(path = atac_frag_path,
                       cells = cells.use)

# Construct the ATAC-seq matrix for the specified cluster
atac.matrix <- FeatureMatrix(
  fragments = list.fragments,
  features = bed.peaks,
  cells = cells.use
)

# Save ATAC-seq matrix
saveRDS(atac.matrix,
	atac_matrix_path)
