## Add Kendall and/or ARC to feature table if needed

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# input files
feature_table = fread(snakemake@input$feature_table_file)
pred_ext = fread(snakemake@input$predictions_extended)
arc_pred_file = snakemake@input$arc_predictions
out_file = snakemake@output$plus_external_features

# check feature table to see if Kendall or ARC.E2G are possible features
input_features = c(feature_table$input_col, feature_table$second_input)
input_features = unique(na.omit(input_features)) # remove NAs, keep unique
to_add = c()
if ("ARC.E2G.Score" %in% input_features) {
	to_add =  c("ARC.E2G.Score", "Kendall")
} else if ("Kendall" %in% input_features) {
	to_add = c("Kendall")
}

# add features if necessary
if (length(to_add)>0) {
  arc_pred = fread(arc_pred_file)
  output <- arc_pred %>%
    select(name, TargetGene, all_of(to_add)) %>%
    left_join(pred_ext, ., by = c("name", "TargetGene"))
} else {
  output = pred_ext
}

# write output
fwrite(output, file = out_file, sep = "\t", na = "NA", quote = FALSE)


