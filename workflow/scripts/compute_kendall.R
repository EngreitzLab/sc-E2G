## Compute Kendall correlation
# Note: Plan to replace this script with an Rcpp version in the future for improved performance

# Load required packages
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(genomation)
  library(foreach)
  library(Signac)
})

## Define functions --------------------------------------------------------------------------------
MyCountDiff = function(y.matrix.sorted){
  n = nrow(y.matrix.sorted)
  cumsum.matrix = apply(y.matrix.sorted, 2, cumsum)
  
  concordant.matrix = cumsum.matrix * (1 - y.matrix.sorted)
  disconcordant.matrix = (1:n - cumsum.matrix) * y.matrix.sorted 
  
  colSums(concordant.matrix - disconcordant.matrix)
}


MyFastKendall = function(x, y.matrix){
  ord = order(x, 
              decreasing = T)
  x.sorted = x[ord]
  y.matrix.sorted = 
    y.matrix[ord, ,drop = F]
  
  n.diff = MyCountDiff(y.matrix.sorted)
  
  x.ties = unique(x.sorted[duplicated(x.sorted)])
  for (x.tie in x.ties) {
    n.diff = 
      n.diff - 
      MyCountDiff(y.matrix.sorted[x.sorted == x.tie, ,drop = F])
  }
  
  l = length(x)
  s = colSums(y.matrix)
  tx = table(x)
  
  n0 = choose(l, 2)
  n1 = sum(choose(tx, 2))
  n2 = (s*(s-1) + (l-s)*(l-s-1))/2
  
  tau_b = n.diff / sqrt((n0 - n1) * (n0 - n2))
  
  return(tau_b)
}

MyFastKendall.E2G <- function(bed.E2G,
                              data.RNA,
                              data.ATAC,
                              colname.gene_name = "gene_name",
                              colname.enhancer_name = "peak_name",
                              colname.output = "Kendall") {

  bed.E2G.filter = 
    bed.E2G[mcols(bed.E2G)[,colname.gene_name] %in% rownames(data.RNA) &
              mcols(bed.E2G)[,colname.enhancer_name] %in% rownames(data.ATAC)] 
  
  data.ATAC = t(BinarizeCounts(data.ATAC))
  
  bed.E2G.output <- foreach(gene.name = unique(mcols(bed.E2G.filter)[,colname.gene_name]),
                     .combine = 'c') %do% {
                       
                       bed.E2G.tmp <- bed.E2G.filter[mcols(bed.E2G.filter)[,colname.gene_name] == gene.name]
                       mcols(bed.E2G.tmp)[, colname.output] = MyFastKendall(as.numeric(data.RNA[gene.name, ]),
                                                                            data.ATAC[, mcols(bed.E2G.tmp)[,colname.enhancer_name], drop = F])
                       bed.E2G.tmp
                     }
  return(bed.E2G.output)
}

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
matrix.atac = read.csv(atac_matix_path,
                       row.names = 1,
                       check.names = F)

# Load scRNA matrix
matrix.rna = read.csv(rna_matix_path,
                      row.names = 1,
                      check.names = F)
matrix.rna = matrix.rna[,colnames(matrix.atac)]

# Compute Kendall correlation
pairs.E2G = MyFastKendall.E2G(pairs.E2G,
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
write.table(df.pairs.E2G,
            file = gzfile(kendall_predictions_path),
            row.names = F,
            quote = F,
            sep = "\t")

