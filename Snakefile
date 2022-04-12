#!/usr/bin/env snakemake

"""
Variance Decomposition pipeline
----------------------

Snakemake pipeline to perform variance decomposition using limix
"""

__version__ = '0.1.0'

container: "docker://henryjt/flexible_variance_decomposition:1.0.0"

rule all:
    input:
        # Variance decomposition
        expand(
           [
            'results/factor_{k}/norm={normalization}__resid={resid}/mtc={mtx_label}/var_decomp.tsv.gz',
            'results/factor_{k}/norm={normalization}__resid={resid}/mtc={mtx_label}/finished_plotting.txt'
           ],
           k=config['PEER_N_FACTORS'],
           normalization=config['VARDEC_NORMALIZATION'],
           resid=config['VARDEC_RESIDUAL_DISTRIBUTION'],
           mtx_label=list(config['VARDEC_MATRICES'].keys())
        ),

        expand(
            'results/factor_{k}/decomposed_variance.tsv.gz',
            k=config['PEER_N_FACTORS']
        )

# File prep ####################################################################
rule prep__moltraits:
    """
    Filter molecular trait tsv to match samples.
    """
    input:
        moltraits='data/moltraits.tsv.gz',
        samples='data/samples.txt'
    output:
        'results/moltraits.tsv.gz'
    params:
        mol_id = config['MOLTRAITS_ID_COL'],
        mol_sample_start = config['MOLTRAITS_SAMPLE_START'],
        script = srcdir('scripts/001-prep_moltraits.py')
    shell:
        '{params.script} '
        '--moltraits {input.moltraits} '
        '--moltrait_id_col {params.mol_id} '
        '--sample_start {params.mol_sample_start} '
        '--samples {input.samples} '
        '--output_file {output}'

rule prep__covariates:
    """
    Prepare covaraite file. If doesn't exist, makes a new covariate file.
    """
    input:
        samples='data/samples.txt'
    output:
        'results/covariates.tsv.gz'
    params:
        cov_file = 'data/covariates.tsv.gz',
        covs = ('--covariates {}'.format(','.join(config['COVARIATES'])) if
            len(config['COVARIATES']) > 0 else ''
        ),
        covs_refs = (
            '--covariate_references {}'.format(
                ','.join(config['COVARIATES_REFERENCES'])
            ) if len(config['COVARIATES_REFERENCES']) > 0 else ''
        ),
        covs_mi = (
            '--covariate_mean_impute {}'.format(
                ','.join(config['COVARIATES_MEAN_IMPUTE'])
            ) if len(config['COVARIATES_MEAN_IMPUTE']) > 0 else ''
        ),
        script = srcdir('scripts/003-prepare_covariates.py')
    shell:
        '{params.script} '
        '--covariate_file {params.cov_file} '
        '--samples {input.samples} '
        '--output_file {output} '
        '{params.covs} '
        '{params.covs_refs} '
        '{params.covs_mi} '

rule prep__covariance_matrices:
    """
    Generate covariances matrices for vectors.
    """
    input:
        samples='data/samples.txt',
        matrix=lambda wildcards: config['MATRICES'][wildcards.matrix]
    output:
        'results/covariance_matrices/{matrix}.tsv.gz'
    params:
        transform=lambda wildcards: (wildcards.matrix in
            config['MATRICES_TRANSFORM_DATAFRAMES']
        ),
        script = srcdir('scripts/005-calculate_covariance_mtx.py')
    run:
        if params.transform:
            shell((
                '{params.script} '
                '--mtx_input {input.matrix} '
                '--samples {input.samples} '
                '--output_file {output} '
            ))
        else:
            shell('cp {input.matrix} {output}')
################################################################################

# PEER Factors #################################################################
rule peer_factors__calculate:
    """
    Calculate PEER factors
    """
    input:
        moltrait = 'results/moltraits.tsv.gz',
        covs = 'results/covariates.tsv.gz'
    output:
        factors = 'results/factor_{k}/peer/moltraits-peer_factors.tsv.gz',
        weights = 'results/factor_{k}/peer/moltraits-peer_weights.tsv.gz',
        precision = 'results/factor_{k}/peer/moltraits-peer_precision.tsv.gz',
        resids = 'results/factor_{k}/peer/moltraits-peer_residuals.tsv.gz',
        resids_inv = 'results/factor_{k}/peer/moltraits-peer_residuals-invnorm.tsv.gz'
    conda:
        'envs/environment-peer.yml'
    params:
        sample_start = config["MOLTRAITS_SAMPLE_START"]+1, # 1-indexed for R
        iterations = config["PEER_ITERATIONS"],
        acct_mean = config["PEER_ACCOUNT_MEAN"],
        inv_norm = config["PEER_INVERSE_NORMALIZE"],
        script = srcdir('scripts/011-calculate_peer_factors.R')
    shell:
        'if [ {wildcards.k} = "NA" ]; then '
        # Make files for downstream analysis.
        '   touch {output.factors}; '
        '   touch {output.weights}; '
        '   touch {output.precision}; '
        '   touch {output.resids}; '
        '   touch {output.resids_inv}; '
        'else '
        # Run PEER
        '   {params.script} '
        '   --phenotypes {input.moltrait} '
        '   --sample_column_start {params.sample_start} '
        '   --hidden_factors {wildcards.k} '
        '   --iterations {params.iterations} '
        '   --covariate_file {input.covs} '
        '   --account_mean {params.acct_mean} '
        '   --inverse_norm {params.inv_norm} '
        '   --output_base results/factor_{wildcards.k}/peer/moltraits; '
        'fi'

rule peer_factors__prep_model_inputs:
    """
    Prepares input for LIMIX models based on PEER factor settings. If
    `use_peer_residuals` == True and `k` != 0, uses residuals for phenotype.
    If not, sym links normal phenotype matrix.

    For covariates, concats PEER factors for downstream analysis.
    """
    input:
        moltraits = 'results/moltraits.tsv.gz',
        covs = 'results/covariates.tsv.gz',
        peer_factors = 'results/factor_{k}/peer/moltraits-peer_factors.tsv.gz',
        peer_resids = 'results/factor_{k}/peer/moltraits-peer_residuals.tsv.gz'
    output:
        moltraits = 'results/factor_{k}/moltraits.tsv.gz',
        covs = 'results/factor_{k}/covariates.tsv.gz'
    params:
        peer_resids = config["PEER_USE_RESIDUALS"],
        script = srcdir('scripts/013-merge_peer_factors_covariates.py')
    run:
        if wildcards.k == 'NA':
            shell((
                'ln -s "$PWD"/{input.moltraits} {output.moltraits}; '
                'ln -s "$PWD"/{input.covs} {output.covs}'
            ))
        else:
            if params.peer_resids:
                # If resid, covariates are included in PEER model. Don't carry
                # downstream. Just copy sample names.
                shell((
                    'ln -s "$PWD"/{input.peer_resids} {output.moltraits}; '
                    # 'ln -s "$PWD"/{input.covs} {output.covs}'
                    'zcat {input.covs} | cut -f1 | gzip > {output.covs}'
                ))
            else:
                # Copy normal pheno, concat factors to covs
                shell((
                    'ln -s "$PWD"/{input.moltraits} {output.moltraits}; '
                    '{params.script} '
                    '--covariates {input.covs} '
                    '--peer_factors {input.peer_factors} '
                    '--output_file {output.covs}'
                ))
################################################################################

# Variance decomposition #######################################################
rule variance_decomposition__decompose:
    """
    Decompose variance using limix.
    """
    input:
        moltraits = 'results/factor_{k}/moltraits.tsv.gz',
        covs = 'results/factor_{k}/covariates.tsv.gz',
        matrices = lambda wildcards: expand(
            'results/covariance_matrices/{matrix}.tsv.gz',
            matrix = config['VARDEC_MATRICES'][wildcards.mtx_label]
        ),
        samples='data/samples.txt'
    output:
        'results/factor_{k}/norm={normalization}__resid={resid}/mtc={mtx_label}/var_decomp.tsv.gz'
    params:
        mol_id = config['MOLTRAITS_ID_COL'],
        mtx_paths = lambda wildcards: ','.join(expand(
            'results/covariance_matrices/{matrix}.tsv.gz',
            matrix = config['VARDEC_MATRICES'][wildcards.mtx_label]
        )),
        mtx_descrip = lambda wildcards: ','.join(
            config['VARDEC_MATRICES'][wildcards.mtx_label]
        ),
        normalization = lambda wildcards: '--normalize {}'.format(
            wildcards.normalization
        ) if wildcards.normalization != 'NA' else '',
        plot_dir = lambda wildcards: 'results/factor_{}/norm={}__resid={}/mtc={}/plots'.format(
            wildcards.k,
            wildcards.normalization,
            wildcards.resid,
            wildcards.mtx_label
        ),
        script = srcdir('scripts/021-variance_decomposition.py')
    shell:
        '{params.script} '
        '--moltraits {input.moltraits} '
        '--moltrait_id_col {params.mol_id} '
        '--matrix_paths {params.mtx_paths} '
        '--matrix_descriptors {params.mtx_descrip} '
        '--additional_covariates {input.covs} '
        '--residual_distribution {wildcards.resid} '
        '--samples {input.samples} '
        '--output_file {output} '
        '--output_plot_dir {params.plot_dir} '
        '{params.normalization}'

rule variance_decomposition__plot:
    """
    Plot variance decomposition
    """
    input:
        'results/factor_{k}/norm={normalization}__resid={resid}/mtc={mtx_label}/var_decomp.tsv.gz'
    output:
        'results/factor_{k}/norm={normalization}__resid={resid}/mtc={mtx_label}/finished_plotting.txt'
    params:
        mtcs = lambda wildcards: ','.join(
            config['VARDEC_MATRICES'][wildcards.mtx_label]
        ),
        out_dir = lambda wildcards: 'results/factor_{}/norm={}__resid={}/mtc={}'.format(
            wildcards.k,
            wildcards.normalization,
            wildcards.resid,
            wildcards.mtx_label
        ),
        script = srcdir('scripts/023-plot_variance_decomposition.py')
    shell:
        '{params.script} '
        '--vardec_file {input} '
        '--effect_labels {params.mtcs} '
        '--output_plot_dir {params.out_dir}; '
        'touch {output}'

rule variance_decomposition__concat:
    """
    Combine the output from the split data
    """
    input:
        expand(
            'results/factor_{{k}}/norm={normalization}__resid={resid}/mtc={mtx_label}/var_decomp.tsv.gz',
            normalization=config['VARDEC_NORMALIZATION'],
            resid=config['VARDEC_RESIDUAL_DISTRIBUTION'],
            mtx_label=config['VARDEC_MATRICES']
        )
    output:
        'results/factor_{k}/decomposed_variance.tsv.gz'
    params:
        files=lambda wildcards: ','.join(expand(
            'results/factor_{k}/norm={normalization}__resid={resid}/mtc={mtx_label}/var_decomp.tsv.gz',
            k=wildcards.k,
            normalization=config['VARDEC_NORMALIZATION'],
            resid=config['VARDEC_RESIDUAL_DISTRIBUTION'],
            mtx_label=config['VARDEC_MATRICES']
        )),
        script = srcdir('scripts/merge_dataframes.py')
    shell:
        '{params.script} '
        '--dataframes {params.files} '
        '--output_file {output}'

################################################################################
