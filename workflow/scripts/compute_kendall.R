## Compute Kendall correlation

# required packages
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(genomation)
  library(foreach)
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

  data.ATAC = t(data.ATAC)
  
  bed.E2G <- foreach(gene.name = unique(mcols(bed.E2G)[,colname.gene_name]),
                     .combine = 'c') %do% {
                       
                       bed.E2G.tmp <- bed.E2G[mcols(bed.E2G)[,colname.gene_name] == gene.name]
                       mcols(bed.E2G.tmp)[, colname.output] = MyFastKendall(as.numeric(data.RNA[gene.name, ]),
                                                                            data.ATAC[, mcols(bed.E2G.tmp)[,colname.enhancer_name], drop = F])
                       bed.E2G.tmp
                     }
  return(bed.E2G)
}



# load candidate E-G pairs
pairs.E2G = readGeneric(snakemake@input$kendall_candidate_e2g_pairs,
                        keep.all.metadata = T,
                        header = T)
# load scRNA matrix
matrix.rna = read.csv(snakemake@input$sc_rna_matix,
                      row.names = 1,
                      check.names = F)

# load scATAC matrix
matrix.atac = read.csv(snakemake@input$sc_atac_matix,
                       row.names = 1,
                       check.names = F)
# compute Kendall correlation
pairs.E2G = MyFastKendall.E2G(pairs.E2G,
                              matrix.rna,
                              matrix.atac,
                              colname.gene_name = "gene_name",
                              colname.enhancer_name = "peak_name",
                              colname.output = "Kendall")

# write output to file
write.table(as.data.frame(pairs.E2G),
          file = gzfile(snakemake@output$kendall_predictions),
	  row.names = F,
	  quote = F,
	  sep = "\t")

