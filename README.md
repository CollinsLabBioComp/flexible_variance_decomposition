Variance Decomposition using LIMIX
==================================

Snakemake pipeline to perform variance decomposition using [LIMIX](https://github.com/limix/limix).

* GitHub repo: https://github.com/CollinsLabBioComp/flexible_variance_decomposition
* Free software: MIT license

# Description

This pipeline is designed to perform variance decomposition in molecular trait (e.g., RNA-seq) data using LIMIX v3.0.4. Briefly, we calculate a covariance matrix from a set of observations (e.g., genetics, differential potential, etc.) and model the matrix as a random effect in LIMIX to estimate the proportion of variance in the molecular trait explained by the effect. As input, it takes molecular trait data and either (1) a data frame of values to convert to a covariance matrix or (2) a pre-calculated covariance matrix.

For detailed description of LIMIX: [LIMIX: genetic analysis of multiple traits](https://www.biorxiv.org/content/10.1101/003905v2)

# Quickstart

Quickstart for deploying this pipeline locally and on a high performance compute cluster.


## 1. Set up the environment

See [environment README](envs/README.md) to set up environment. Once the environment is set up, activate the `conda` environment:

```bash
source activate limix-vardec
```

Alternatively, if using singularity or docker, one can pull the image from [henryjt/flexible_variance_decomposition:1.0.0](https://hub.docker.com/layers/202247912/henryjt/flexible_variance_decomposition/1.0.0/images/sha256-9fa1a061f8c5f2f9cec93a08ac61a20652e068b41924a4cdbaaf6ec1349203c8?context=repo).

## 2. Prepare the input files

Generate and/or edit input files for the pipeline.

As input, the pipeline expects the following files in the `data/` directory in the location you are running the pipeline:
1. `moltraits.tsv.gz`: A TSV file containing molecular trait data where samples are columns and rows are features. All metadata information should be before the
    samples start. For example:
    | chr | start | end | gene | sample_1 | sample_2 | ... | sample_n |
    | --- | ----- | --- | ---- | -------- | -------- | --- | -------- |
    | chr11 | 2159779 | 2161221 | INS | 40.5 | 241.5 | ... | 591.1 |
    | chr2 | 162142882 | 162152404 | GCG | 72.5 | 10.5 | ... | 1000.1 |
2. `samples.txt`: Text file containing all samples to use for analysis. Each line should be a sample ID and should be contained within `moltraits.tsv.gz` and `covariates.tsv.gz`
3. `covariates.tsv.gz`: TSV file containing covariates to include in the model. First column should correspond to sample IDs.

Examples of these files can be found in `demo/`.


## 3. Run pipeline

**NOTE**: All input file paths should be full paths.

To run:
```bash
snakemake \
   --snakefile "/path/to/repo/dir/Snakefile" \
   --configfile "/path/to/config/config_analysis.json"
```

Examples:
* [Local with docker](demo/run_variance_decomposition__docker.sh)
* [Local with conda](demo/run_variance_decomposition__local.sh)
* [Cluster with conda](demo/run_variance_decomposition__sge.sh)

# Notes
* The primary LIMIX documentation is no longer being supported. To find more information, see the [temporary documentation](https://limix-tempdoc.readthedocs.io/en/latest/index.html).


Authors: Henry Taylor
