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
     --latency-wait 600 \
     --jobs 999 \
     --cluster-config ${SNK_DIR}/lib/configs/config_cluster_sge.json \
     --cluster ${SNK_DIR}/lib/wrappers/cluster/sge.py \
     $1
