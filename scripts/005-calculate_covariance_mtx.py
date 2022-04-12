#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import scipy as sp
import numpy as np
import pandas as pd


def calculate_sample_covariance(mtx):
    """
    Calculated a sample covariance matrix (X)(X.T).

    Input:
        * mtx: matrix of data. Rows are samples, columns are variables.

    Output:
        * Unnormalized matrix
    """
    # this is what that LIMIX function actually does
    mtx -= mtx.mean(0)    # subtract the mean for each feature across samples
    mtx /= mtx.std(0)     # divide by std dev for each feature across samples
    cov_mtx = sp.dot(mtx, mtx.T)

    # Normalize data for downstream analyses
    # NOTE: Don't normalize here -- normalize after we combine across chroms.
    # norm_cov_mtx = limix.qc.normalise_covariance(cov_mtx)

    return cov_mtx


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="Calculate covariance matrix from a dataframe."
    )

    parser.add_argument(
        '-mtx', '--mtx_input',
        action='store',
        dest='mtx',
        default='',
        help='''
        TSV containing values to create covariance matrix. First column should
        be sample IDs and the rest columns should be input into covariance
        matrix. For example:
        sample_ids
        variable_1
        variable_2
        ...
        variable_n
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
        default='sample_covariance.tsv.gz',
        help='Output file.'
    )
    options = parser.parse_args()

    # Load in data
    df = pd.read_csv(options.mtx, sep='\t', header=0, index_col=0)
    samples = [
        str(x) for x in # in case the samples are numbers. Pandas will read headers as string
        pd.read_csv(options.samples, sep='\n', header=None)[0].to_list()
    ]

    # First check to make sure all samples are found in TSV
    missing_samp = np.setdiff1d(samples, df.index.to_list()).tolist()
    if len(missing_samp) > 0:
        raise ValueError(
            'Some samples specified in {} are not in {}: {}'.format(
                options.samples,
                options.mtx,
                ', '.join(missing_samp)
            )
        )

    var_mtx = np.matrix(df.loc[samples, ])
    smpl_cov_mtx = calculate_sample_covariance(var_mtx)

    # Write to TSVs
    pd.DataFrame(smpl_cov_mtx).to_csv(
        options.of,
        index=False,
        header=False,
        sep='\t',
        compression='gzip'
    )


if __name__ == '__main__':
    main()
