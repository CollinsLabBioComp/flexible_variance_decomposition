#!/bin/sh

## Activate QTL conda environment
source activate limix-vardec

## Set any necessary env variables
export SNK_DIR=$(realpath ..)

## ADD ANY SPECIFIC ENV NEEDS

snakemake \
     --snakefile ${SNK_DIR}/Snakefile \
     --configfile config_analysis.json \
     --use-conda \
     --printshellcmds \
     $1
