# Load required packages
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(genomation)
})

## Define functions --------------------------------------------------------------------------------

# Add Kendall correlation as a new colunm of ABC prediction
AddMax = function(bed.E2G.A,
                  bed.E2G.B,
                  colname.in,
                  colname.out){
  overlaps.res = findOverlaps(bed.E2G.A,
                              bed.E2G.B)
  df.overlaps = as.data.frame(overlaps.res)
  rm(overlaps.res)
  df.overlaps.filter = 
    df.overlaps[bed.E2G.A$TargetGene[df.overlaps$queryHits] == 
                  bed.E2G.B$TargetGene[df.overlaps$subjectHits],]
  rm(df.overlaps)
  
  df.overlaps.filter[,colname.in] = mcols(bed.E2G.B)[df.overlaps.filter$subjectHits,colname.in]
  df.overlaps.filter.sort = df.overlaps.filter[order(df.overlaps.filter[,colname.in],decreasing = T),]
  df.overlaps.filter.sort2 = df.overlaps.filter.sort[order(df.overlaps.filter.sort$queryHits),]
  df.overlaps.filter.sort2.unique = df.overlaps.filter.sort2[!duplicated(df.overlaps.filter.sort2$queryHits),]
  
  
  
  mcols(bed.E2G.A)[df.overlaps.filter.sort2.unique$queryHits,
                   colname.out] = 
    mcols(bed.E2G.B)[df.overlaps.filter.sort2.unique$subjectHits,
                     colname.in]
  bed.E2G.A
}

# Integrate ABC and Kendall into ARC-E2G
IntergrateABC = function(bed.E2G,
                         colname.ABC,
                         colname.feature,
                         colname.output){
  
  index.exclude_na = !is.na(mcols(bed.E2G)[,colname.ABC]) &
    mcols(bed.E2G)[,colname.ABC] > 0 &
    !is.na(mcols(bed.E2G)[,colname.feature])
  
  bed.E2G.filter = bed.E2G[index.exclude_na]
  
  sd.lnABC = sd(log(mcols(bed.E2G.filter)[,colname.ABC]))
  sd.feature = sd(mcols(bed.E2G.filter)[,colname.feature])
  
  lm.res = 
    lm(scale(log(mcols(bed.E2G.filter)[,colname.ABC])) ~
         scale(mcols(bed.E2G.filter)[,colname.feature]))
  
  beta = 1/coef(lm.res)[2]
  
  mcols(bed.E2G)[,colname.output] = mcols(bed.E2G)[,colname.ABC]
  mcols(bed.E2G)[index.exclude_na,colname.output] =
    exp(log(mcols(bed.E2G)[index.exclude_na,colname.ABC]) +
          mcols(bed.E2G)[index.exclude_na,colname.feature] * sd.lnABC / sd.feature / beta)
  
  return(bed.E2G)
}

# Import parameters from Snakemake
abc_predictions_path = snakemake@input$abc_predictions
kendall_predictions_path = snakemake@input$kendall_predictions
arc_predictions_path = snakemake@output$arc_predictions

# Load ABC data
pairs.E2G.ABC = readGeneric(abc_predictions_path,
                            keep.all.metadata = T,
                            header = T)
# Load Kendall data
pairs.E2G.Kendall = readGeneric(kendall_predictions_path,
                                keep.all.metadata = T,
                                header = T)
# Add Kendall correlation as a new column of ABC prediction
pairs.E2G.ABC = AddMax(pairs.E2G.ABC,
                       pairs.E2G.Kendall,
                       "Kendall",
                       "Kendall")

# Compute ARE-E2G
pairs.E2G.ABC = IntergrateABC(pairs.E2G.ABC,
                              "ABC.Score",
                              "Kendall",
                              "ARC.E2G.Score")

# Write output to file
write.table(as.data.frame(pairs.E2G.ABC),
            file = gzfile(arc_predictions_path),
            row.names = F,
            quote = F,
            sep = "\t")

