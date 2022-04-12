#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import pandas as pd
import scipy as sp
import limix.vardec as lm
from pathlib import Path
from math import ceil
import plotnine as plt9
from optimix import OptimixError

# Import useful functions from other scripts
import _normalize_phenotype as norm_pheno


def decompose_variance(
    phenotype,
    residual_distribution,
    covariance_mtx_dict,
    covs
):
    var_decomp = lm.VarDec(
        y=phenotype,
        lik=residual_distribution,
        M=covs
    )
    for key in covariance_mtx_dict.keys():
        var_decomp.append(covariance_mtx_dict[key], key)
    var_decomp.append_iid("noise")
    var_decomp.fit(verbose=False)
    return var_decomp


def plot_variance_decomposition(
    data,
    labels
):
    df = pd.DataFrame({
        'values': data,
        'labels': pd.Categorical(labels, categories=labels, ordered=True)
    })
    plt = (
        plt9.ggplot(df, plt9.aes(x='labels', y='data', fill='labels'))
        + plt9.geom_col()
        + plt9.labs(
            x='Effect',
            y='Variance Explained',
            title='Variance Decomposition'
        )
        + plt9.theme_bw()
        + plt9.theme(legend_position='none')
    )
    return plt


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="Variance decomposition."
    )

    parser.add_argument(
        '-mt', '--moltraits',
        action='store',
        dest='mt',
        required=True,
        help='''
        TSV file containing moltrait information. Moltrait ID column should be
        specificed with `moltrait_id_col`. Sample information should start at
        column `sample_start`. All moltrait information should be before the
        samples start. For example:
        chr
        start
        end
        <moltrait_info_1>
        <phenotype-id>: specified with `moltrait_id_col` flag
        <sample_names>: start of sample information
        '''
    )

    parser.add_argument(
        '-mid', '--moltrait_id_col',
        action='store',
        dest='mid',
        required=True,
        help='''
        Column name specifying moltrait ID
        '''
    )

    parser.add_argument(
        '-mtx_paths', '--matrix_paths',
        action='store',
        dest='mtx_paths',
        default=None,
        help='''
        Comma-separated list of additional matrices to add in variance
        decomposition analyses.
        '''
    )

    parser.add_argument(
        '-mtx_descr', '--matrix_descriptors',
        action='store',
        dest='mtx_descr',
        default=None,
        help='''
        Comma-separated list of descriptions of `mtx_paths`.
        '''
    )

    parser.add_argument(
        '-covs', '--additional_covariates',
        action='store',
        dest='covs',
        default=None,
        help='''
        Any additional covariates to use. Rows shold be samples. Columns should
        be numeric covariates to include. First column should be sample IDs.
        *Note*: All columns will be included as covariates.
        '''
    )

    parser.add_argument(
        '-rd', '--residual_distribution',
        action='store',
        dest='rd',
        default='normal',
        help='''
        Expected distribution of residual effect sizes in the fitted linear
        mixed model. Can be one of the following:
        "normal", "bernoulli", "probit", "binomial", "poisson"
        (default: %(default)s).
        '''
    )

    parser.add_argument(
        '-n', '--normalize',
        action='store',
        dest='normalize',
        default=None,
        help='''
        Normalize sample expression to a give distribution PER gene. Options:
        'sklearn_box_cox', 'sklearn_rank_normal' or 'rank_normal'
        (default: %(default)s).
        '''
    )

    parser.add_argument(
        '-s', '--samples',
        action='store',
        dest='samples',
        required=True,
        help='''
        Text file containing samples to retain. Each line should be a sample
        name.
        '''
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='decomposed_variance.tsv.gz',
        help='Output file for all traits.'
    )

    parser.add_argument(
        '-opd', '--output_plot_dir',
        action='store',
        dest='opd',
        default='plots',
        help='Output directory for plots'
    )
    options = parser.parse_args()

    # Grab data
    moltraits = pd.read_csv(options.mt, sep='\t', header=0)
    samples = [
        str(x) for x in # in case the samples are numbers. Pandas will read headers as string
        pd.read_csv(options.samples, sep='\n', header=None)[0].to_list()
    ]

    # Get additional matrices
    mtcs = {}
    mtx_paths = options.mtx_paths.split(',')
    mtx_descripts = options.mtx_descr.split(',')
    for i,x in enumerate(mtx_paths):
        mtcs[mtx_descripts[i]] = np.asmatrix(
            pd.read_csv(x, sep='\t', header=None, index_col=None)
        )
        print(mtcs[mtx_descripts[i]])


    covs = pd.read_csv(options.covs, sep='\t', header=0, index_col=0)
    if covs.size == 0: # check to see if no data in covariates
        covs = None

    # Make sure plot path exists
    # Also make sure directory for var decomp plots exists
    Path(options.opd).mkdir(parents=True, exist_ok=True)

    df_columns = (
        ['moltrait_id', 'normalization', 'method', 'effects_analyzed'] +
        list(mtcs.keys()) +
        ['noise', 'mean_expr']
    )
    vardec_df = pd.DataFrame([], columns=df_columns)

    for trait in moltraits.iterrows():
        # iterrows returns (index, data) -- grab data
        phenotype = trait[1]
        pheno_id = phenotype[options.mid]
        print('Decomposing variance for trait: {}...'.format(pheno_id))

        # Get phenotype matrix
        pheno_mtx = phenotype[samples]  # grab samples
        mean_expr = np.mean(pheno_mtx.values)  # Get mean expr before norm
        pheno_mtx = norm_pheno.correct_phenotype(
            phenotype_df=pheno_mtx,
            normalize_method=options.normalize
        )

        # Decompose variance
        try:
            # LIMIX normalizes the covariance matrices, so don't normalize
            # before call
            var_decomp = decompose_variance(
                phenotype=pheno_mtx,
                residual_distribution=options.rd,
                covariance_mtx_dict=mtcs,
                covs=covs
            )
            print(var_decomp)
            print(var_decomp.covariance)
            print(var_decomp.covariance[0])
            print(var_decomp.covariance[1])

            pheno_vardec_df = pd.DataFrame(
                [[
                    pheno_id, options.normalize, options.rd,
                    ','.join(list(mtcs.keys())), mean_expr
                ]],
                columns=[
                    'moltrait_id', 'normalization', 'method',
                    'effects_analyzed', 'mean_expr'
                ]
            )
            for i,x in enumerate(mtcs.keys()):
                pheno_vardec_df[x] = var_decomp.covariance[i].scale
            pheno_vardec_df['noise'] = var_decomp.covariance[len(mtcs)].scale

            # Plot result
            plt_vals = list(mtcs.keys()) + ['noise']
            plt = plot_variance_decomposition(
                pheno_vardec_df[plt_vals].values[0],
                plt_vals
            )
            plt.save(
                filename='{}/{}_var_decomp.png'.format(options.opd, pheno_id),
                dpi=320
            )

        except OptimixError:
            print((
                'Limix variance decomposition model could not converge. '
                'Returning all `NA`s'
            ))
            empty_val = [
                pheno_id, options.normalize, options.rd,
                ','.join(list(mtcs.keys()))
            ]
            empty_val = empty_val + (['NA'] * (len(mtcs)+1)) + [mean_expr]
            pheno_vardec_df = pd.DataFrame([empty_val], columns=df_columns)

        except ValueError:
            print((
                'Limix variance decomposition model could not converge. '
                'Returning all `NA`s'
            ))
            empty_val = [
                pheno_id, options.normalize, options.rd,
                ','.join(list(mtcs.keys()))
            ]
            empty_val = empty_val + (['NA'] * (len(mtcs)+1)) + [mean_expr]
            pheno_vardec_df = pd.DataFrame([empty_val], columns=df_columns)

        vardec_df = vardec_df.append(pheno_vardec_df, ignore_index=True)
        vardec_df = vardec_df.replace(np.nan,'NA')
        print('Done.')

    vardec_df.to_csv(
        options.of,
        index=False,
        sep='\t',
        compression='gzip',
    )


if __name__ == '__main__':
    main()
