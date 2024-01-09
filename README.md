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
conda create -n mamba -c conda-forge mamba python=3.11
conda activate mamba
mamba install -c conda-forge -c bioconda snakemake=7
```

## Apply model
Before running this workflow, users should perform clustering to define cell clusters through regular single-cell analysis (such as Seurat & Signac).
Required input data includes (refer to the example data in the `example_chr22_multiome` folder):
1. Combined meta_data, including columns for cell_name, sample, cluster, and barcode.
2. Combined RNA matrix, normalized by the total UMI count of each cell. It can be with or without log transformation.
3. Indexed fragment files for each sample.

To configure the pipeline:
- Modify `config/config.yaml` to specify paths for meta_data, rna_matrix, and results_dir.
- Modify `config/config_multiome_samples.tsv` to specify the fragment file path for each sample.
- Modify `config/config_cell_clusters.tsv` to specify the Hi-C file path, Hi-C data type, Hi-C resolution, TSS coordinates, and gene coordinates for each cell cluster.



Running the pipeline:
```
snakemake -j1 --use-conda
```
This command make take a while the first time you run it, as it needs to build the conda environments. 
But if it takes more than 1 hour, that's usually a bad sign, meaning that you're not using mamba.

Output will show up in the `results/` directory by default

## Train model

**Important: Only train models for biosamples matching the corresponding CRISPR data (in this case, K562)**

Modify `config/config_training.yaml` with your model and dataset configs
- `model_config` has columns: model, dataset, ABC_directory, feature_table, polynomial (do you want to use polynomial features?) 
Note that trained models generated using polynomial features cannot directly be used in the **Apply model** workflow
- `dataset_config` has rows representing each "dataset"  in `model_config`, where each "dataset" must correspond to a "biosample" in `dataset_config`
    - Each "dataset"/"biosample" **must** have an ATAC_matrix, RNA_matrix, and candidate_pairs
    - If an ABC_directory is not specified for a dataset, its entry in `dataset_config` must also contain the required ABC biosample parameters
    - TO DO: specify how to generate and formats for Kendall parameters

Running the pipeline:
```
snakemake -s workflow/Snakefile_training -j1 --use-conda
```