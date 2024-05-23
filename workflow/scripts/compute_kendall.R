## Compute Kendall correlation

# Load required packages
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(genomation)
  library(foreach)
  library(Signac)
  library(Seurat)
  library(Rcpp)
  library(data.table)
  library(Matrix)
})

## Define functions --------------------------------------------------------------------------------

# Calculate the difference between concordant and disconcordant pairs from a sorted logical matrix
cppFunction('
NumericVector count_diff(LogicalMatrix y_matrix_sorted) {
    int n = y_matrix_sorted.nrow();
    int m = y_matrix_sorted.ncol();
    NumericVector result(m);
    for (int j = 0; j < m; j++) {
        long long concordant = 0;
        long long disconcordant = 0;
        long long cumsum = 0;
        for (int i = 0; i < n; i++) {
            bool tmp = y_matrix_sorted(i, j);
            cumsum += tmp;
            if (tmp) {
                disconcordant += (i + 1 - cumsum);
            } else {
                concordant += cumsum;
            }
        }
        result[j] = static_cast<double>(concordant - disconcordant);
    }
    return result;
}
')

# Compute Kendall correlation between a single gene and multiple enhancers
kendall_one_gene = function(x, y.matrix){
  
  # Sort x in decreasing order and accordingly sort y.matrix
  ord = order(x, 
              decreasing = T)
  x.sorted = x[ord]
  y.matrix.sorted = 
    y.matrix[ord, ,drop = F]
  
  # Calculate initial differences between concordant and disconcordant pairs
  n.diff = count_diff(as.matrix(y.matrix.sorted))
  
  # Adjust differences for ties in x
  x.ties = unique(x.sorted[duplicated(x.sorted)])
  for (x.tie in x.ties) {
    n.diff = 
      n.diff - 
      count_diff(as.matrix(y.matrix.sorted[x.sorted == x.tie, ,drop = F]))
  }
  
  # Calculate Kendall's tau-b coefficient
  l = length(x)
  s = colSums(y.matrix)
  tx = table(x)
  
  n0 = choose(l, 2)
  n1 = sum(choose(tx, 2))
  n2 = (s*(s-1) + (l-s)*(l-s-1))/2
  
  tau_b = n.diff / sqrt((n0 - n1) * (n0 - n2))
  
  return(tau_b)
}


# Compute Kendall correlation between a mutliple genes and multiple enhancers
kendall_mutliple_genes = function(bed.E2G,
                                  data.RNA,
                                  data.ATAC,
                                  colname.gene_name = "gene_name",
                                  colname.enhancer_name = "peak_name",
                                  colname.output = "Kendall") {
  
  # Filter E2G pairs based on presence in RNA and ATAC data
  bed.E2G.filter = 
    bed.E2G[mcols(bed.E2G)[,colname.gene_name] %in% rownames(data.RNA) &
              mcols(bed.E2G)[,colname.enhancer_name] %in% rownames(data.ATAC)] 
  

  
  # Compute Kendall correlation for each gene
  bed.E2G.output <- foreach(gene.name = unique(mcols(bed.E2G.filter)[,colname.gene_name]),
                            .combine = 'c') %do% {
                              
                              bed.E2G.tmp <- bed.E2G.filter[mcols(bed.E2G.filter)[,colname.gene_name] == gene.name]
                              
                              mcols(bed.E2G.tmp)[, colname.output] = 
                                kendall_one_gene(as.numeric(data.RNA[gene.name, ]),
                                                 t(data.ATAC[mcols(bed.E2G.tmp)[,colname.enhancer_name], , drop = F]))
                              bed.E2G.tmp
                            }
  return(bed.E2G.output)
}
## -------------------------------------------------------------------------------------------------

# Import parameters from Snakemake
kendall_pairs_path = snakemake@input$kendall_pairs_path
atac_matix_path = snakemake@input$atac_matix
rna_matix_path = snakemake@input$rna_matix
kendall_predictions_path = snakemake@output$kendall_predictions

# Load candidate E-G pairs
pairs.E2G = readGeneric(kendall_pairs_path,
                        keep.all.metadata = T,
                        header = T)

# Load scATAC matrix
matrix.atac_count = readRDS(atac_matix_path)
matrix.atac = BinarizeCounts(matrix.atac_count)
rm(matrix.atac_count)

# Load scRNA matrix
matrix.rna_count = read.csv(rna_matix_path,
                      row.names = 1,
                      check.names = F)
matrix.rna_count = Matrix(as.matrix(matrix.rna_count), sparse = TRUE)
matrix.rna_count = matrix.rna_count[,colnames(matrix.atac)]

# Normalize scRNA matrix
matrix.rna = NormalizeData(matrix.rna_count)
rm(matrix.rna_count)

# Compute Kendall correlation
pairs.E2G = kendall_mutliple_genes(pairs.E2G,
                                   matrix.rna,
                                   matrix.atac,
                                   colname.gene_name = "TargetGene",
                                   colname.enhancer_name = "PeakName",
                                   colname.output = "Kendall")

# Write output to file
df.pairs.E2G = 
  as.data.frame(pairs.E2G)[,c("seqnames",
                              "start",
                              "end",
                              "TargetGene",
                              "PeakName",
                              "PairName",
                              "Kendall")]
colnames(df.pairs.E2G) = 
  c("chr",
    "start",
    "end",
    "TargetGene",
    "PeakName",
    "PairName",
    "Kendall")
fwrite(df.pairs.E2G,
       file = kendall_predictions_path,
       row.names = F,
       quote = F,
       sep = "\t")