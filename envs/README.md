
# Setting up environment

We use [conda](https://docs.conda.io/en/latest/) for managing necessary packages. **Note**: we recommend using [mamba](https://mamba.readthedocs.io/en/latest/) instead of base conda to speed up download times. If mamba is not available, conda should still work.

Alternatively, we also provide a pre-compiled Docker image to use with [Snakemake's singularity function](https://snakemake.readthedocs.io/en/v5.27.4/snakefiles/deployment.html#running-jobs-in-containers).

## Conda

Install the required packages via conda/mamba:
```bash
# The repo directory.
REPO_MODULE="${HOME}/repo/path/to/this/pipeline"

# Install environment using Mamba. Replace 'mamba' with 'conda' if mamba not available.
mamba env create --name limix-vardec --file ${REPO_MODULE}/envs/environment.yml

# Activate the new Conda environment.
source activate limix-vardec

# To update environment file:
#conda env export --no-builds | grep -v prefix | grep -v name > environment.yml
```

## Docker

Alternatively, we have developed a Docker image using [Dockerfile](Dockerfile). To use the [pre-generated docker image](https://hub.docker.com/layers/202247912/henryjt/flexible_variance_decomposition/1.0.0/images/sha256-9fa1a061f8c5f2f9cec93a08ac61a20652e068b41924a4cdbaaf6ec1349203c8?context=repo), use snakemake's singularity integration: [demo file](../demo/run_variance_decomposition__docker.sh).

**Note**: We include Snakemake in the conda environment. If using the Docker setup, the user will need to have an independent installation of Snakemake (v5.27.4).
