#!/bin/sh

## Set any necessary env variables
export SNK_DIR=$(realpath ..)

## ADD ANY SPECIFIC ENV NEEDS

/path/to/snakemake/bin/snakemake \
     --snakefile ${SNK_DIR}/Snakefile \
     --configfile config_analysis.json \
     --use-singularity \
     --singularity-args "-B ${SNK_DIR}:${SNK_DIR}" \
     --printshellcmds \
     --latency-wait 600 \
     --jobs 999 \
     --cluster-config ${SNK_DIR}/lib/configs/config_cluster_sge.json \
     --cluster ${SNK_DIR}/lib/wrappers/cluster/sge.py \
     $1
