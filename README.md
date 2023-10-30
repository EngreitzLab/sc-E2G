# ENCODE E2G Features
Generate [Encode E2G](https://github.com/karbalayghareh/ENCODE-E2G) input features based on [ABC](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction) output

Supports running ENCODE E2G with a pre-trained model. Put the pretrained model in the config, under `models_dir`

## Usage

Modify `config/config.yml` with your ABC biosamples config

```
mamba env create -f workflow/envs/encode_e2g_features.yml
conda activate encode_e2g_features
snakemake -j1
```

Output will show up in the `results/` directory
