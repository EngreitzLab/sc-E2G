# sc-E2G
Pipeline for running single cell ENCODE E2G.

Input: Multiome dataset (ATAC and RNA)
Output: Enhancer -> Gene Predictions

The pipeline consists of the following components:
1. Compute ABC model predictions
2. Compute Kendall Correlation
3. Based on 1 & 2, construct a feature file to be used as input to the ML model
4. Have ML model make predictions (assigning a score to each E-G pair)

## Usage

Modify `config/config_biosamples.tsv` with your multiome data

```
mamba env create -f workflow/envs/sc_e2g.yml
conda activate sc_e2g
snakemake -j1
```

Output will show up in the `results/` directory
