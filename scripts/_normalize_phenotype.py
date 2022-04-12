#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
import plotnine as plt9
import scipy.stats as ss
import sklearn


def _plot_phenotype(phenotype, out_base):
    df = pd.DataFrame(
        data={'sample': phenotype}
    )
    # Plot distribution
    plt = (
        plt9.ggplot(df, plt9.aes(x='sample'))
        + plt9.geoms.geom_histogram(bins=8)
        + plt9.labs(x='Trait expression', y='Frequency')
    )
    plt.save(
        filename='{}_histogram.png'.format(out_base),
        dpi=320
    )

    # First plot QQ plot
    plt = (
        plt9.ggplot(df, plt9.aes(sample='sample'))
        + plt9.stat_qq()
        + plt9.stat_qq_line()
    )
    plt.save(
        filename='{}_qqplot.png'.format(out_base),
        dpi=320
    )
    return phenotype


def _rank_to_normal(rank, n, c=3.0/8):
    # Standard quantile function
    x = (rank - c) / (n - 2*c + 1)
    return ss.norm.ppf(x)


def _normalize_expression(phenotype, normalization_methed='rank_normal'):
    # need to cast to point -- currently an array of objects
    if normalization_methed == 'sklearn_box_cox':
        norm_mtx = sklearn.preprocessing.power_transform(
            np.array(phenotype).reshape(-1, 1),
            method='box-cox',
            standardize=True,
            copy=True
        )[:0]
    elif normalization_methed == 'sklearn_rank_normal':
        norm_mtx = sklearn.preprocessing.quantile_transform(
            np.array(phenotype).reshape(-1, 1),
            axis=0,
            output_distribution='normal',
            random_state=126,
            copy=True
        )[:0]
    elif normalization_methed == 'rank_normal':
        np.random.seed(1270234881)  # set seed
        # Shuffle by index
        orig_idx = phenotype.index
        alg_input = phenotype.loc[np.random.permutation(orig_idx)].copy()

        # Get rank, ties are determined by their position in the series (hence
        # why we randomised the series)
        rank = ss.rankdata(alg_input, method='ordinal')
        rank = pd.Series(rank, index=alg_input.index)

        # Convert rank to normal distribution
        norm_mtx = rank.apply(_rank_to_normal, n=len(rank), c=3.0/8)[orig_idx]
    else:
        print((
            'Normalization method is not `box_cox`, `sklearn_rank_normal`, '
            'or `rank_normal`. Returning the unnormalized data.'
        ))
    return norm_mtx


def _zscore_standardize(phenotype_df):
    '''
    Standardize the phenotype using Z-score standardization.

    Params:

    phenotype_df: Pandas Series for a specific phenotype.
    '''
    standard_mtx = phenotype_df.copy()
    standard_mtx -= standard_mtx.mean(0)
    standard_mtx /= standard_mtx.std(0)
    return standard_mtx


def correct_phenotype(phenotype_df, normalize_method=None):
    '''
    Normalizes a phenotype using the method defined by `normalize_method`. If
    `normalize_method` is None, standardize the phenotype.

    Params

    phenotype_df: Pandas Series for a specific phenotype.

    normalize_method: str. method to normalize with. Options: `box_cox`,
    `sklearn_rank_normal`, `rank_normal`, or None.
    '''
    if normalize_method is not None:
        pheno_mtx = _normalize_expression(phenotype_df, normalize_method)
    else:
        pheno_mtx = _zscore_standardize(phenotype_df)
    return pheno_mtx
