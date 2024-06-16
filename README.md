CircleCI [![CircleCI](https://dl.circleci.com/status-badge/img/gh/EngreitzLab/sc-E2G.svg?style=svg)](https://app.circleci.com/pipelines/github/EngreitzLab/sc-E2G?branch=main)

# sc-E2G
Pipeline for running single cell ENCODE E2G.

Input: Multiome dataset (ATAC and RNA)
Output: Enhancer -> Gene Predictions (for each cell cluster)

The pipeline consists of the following components:
1. Compute ABC model predictions for each cell cluster
2. Generate E2G features from ABC predictions
3. Compute Kendall correlation and/or ARC-E2G score for each cell cluster
4. Combine 2 & 3 to construct a feature file to be used as input to predictive model
5. (optional) Train predictive model using CRISPR-validated E-G pairs from K562 dataset
6. Apply trained model to make predictions by assigning a score to each E-G pair

## Set up
Clone the repo and set it up for submdoule usage
```
git clone --recurse-submodules https://github.com/EngreitzLab/sc-E2G.git
git config --global submodule.recurse true
```

When running for the first time, the conda environments have to be setup.
For speed, it's recommended that your current environment has mamba installed

```
conda config --set channel_priority flexible  # Make sure your conda config uses flexible channel packaging to prevent unsatisfiable errors
conda create -n mamba -c conda-forge mamba "python<=3.11"
conda activate mamba
mamba install -c conda-forge -c bioconda snakemake=7
```

## Apply model
Before running this workflow, users should perform clustering to define cell clusters through regular single-cell analysis (such as Seurat & Signac).
Required input data includes (refer to the example data in the `resources/example_chr22_multiome_cluster` folder):
1. RNA count matrix (gene x cell) for each cell cluster.
2. Fragment files and their corresponding *.tbi index files for each cell cluster.

Note:
- The RNA count matrix should be in either .csv.gz or .h5ad format.
- The 4th column of the fragment file should correspond to the cell names in the RNA count matrix.
- The fragment file should be sorted by coordinates. If it is not sorted, you can use the sortBed tool from bedtools: `sortBed atac_fragments.unsorted.tsv > atac_fragments.tsv`.
- To create a .tbi index, use bgzip from HTSlib (https://github.com/samtools/htslib) instead of gzip to compress the fragment file: `bgzip atac_fragments.tsv`, then generate the corresponding .tbi index file using `tabix -p bed atac_fragments.tsv.gz`.

To configure the pipeline:
- Modify `config/config.yaml` to specify paths for results_dir.
- Modify `config/config_cell_clusters.tsv` to specify the RNA matrix path, fragment file path, Hi-C file path, Hi-C data type, Hi-C resolution, TSS coordinates, and gene coordinates for each cell cluster.



Running the pipeline:
```
snakemake -j1 --use-conda
```
This command make take a while the first time you run it, as it needs to build the conda environments. 
But if it takes more than 1 hour, that's usually a bad sign, meaning that you're not using mamba.

Output will show up in the `results/` directory by default, with the structure `results/cell_cluster/model_name/encode_e2g_predictions.tsv.gz`. The score column to use is `ENCODE-rE2G.Score.qnorm`.

## Train model

**Important: Only train models for biosamples matching the corresponding CRISPR data (in this case, K562)**

Modify `config/config_training.yaml` with your model and cell_cluster configs
- `model_config` has columns: model, dataset, ABC_directory, feature_table, polynomial (do you want to use polynomial features?) 
Note that trained models generated using polynomial features cannot directly be used in the **Apply model** workflow
- `cell_cluster_config` has rows representing each "dataset"  in `model_config`, where each "dataset" must correspond to a "cluster" in `cell_cluster_config`
    - Each "dataset"/"cluster" **must** have an RNA matrix, fragment file
    - If an ABC_directory is not specified for a dataset, its entry in `cell_cluster_config` must also contain the required ABC biosample parameters
    - TO DO: specify how to generate and formats for Kendall parameters
- To apply a trained model, it must contain the following files: 1) `model.pkl`, 2) `feature_table.tsv`, 3) `threshold_.XX` , 4) `qnorm_reference.tsv.gz` (single column with header `ENCODE-rE2G.Score` that contains raw scores for genomewide predictions)

Running the pipeline:
```
snakemake -s workflow/Snakefile_training -j1 --use-conda
```
Output will show up in the `results_training/` directory by default
