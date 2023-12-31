# sc-E2G
Pipeline for running single cell ENCODE E2G.

Input: Multiome dataset (ATAC and RNA)
Output: Enhancer -> Gene Predictions

The pipeline consists of the following components:
1. Compute ABC model predictions
2. Generate E2G features from ABC predictions
3. Compute Kendall Correlation
5. Combine 2 & 3 to construct a feature file to be used as input to the ML model
6. Have ML model make predictions (assigning a score to each E-G pair)

## Usage
Clone the repo and set it up for submdoule usage
```
git clone --recurse-submodules https://github.com/EngreitzLab/sc-E2G.git
git config --global submodule.recurse true
```

Modify `config/config_biosamples.tsv` with your multiome data

When running for the first time, the conda environments have to be setup.
For speed, it's recommended that your current environment has mamba installed

```
conda config --set channel_priority flexible  # Make sure your conda config uses flexible channel packaging to prevent unsatisfiable errors
conda create -n mamba -c conda-forge mamba
conda activate mamba
mamba install -c conda-forge -c bioconda snakemake
```

Running the pipeline:
```
snakemake -j1 --use-conda
```
This command make take a while the first time you run it, as it needs to build the conda environments. 
But if it takes more than 1 hour, that's usually a bad sign, meaning that you're not using mamba.

Output will show up in the `results/` directory by default

## Training the Model

(Maya to fill out)
