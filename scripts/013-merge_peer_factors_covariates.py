#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import pandas as pd
import sys


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="Merge PEER factors into covariates file."
    )

    parser.add_argument(
        '-covs', '--covariates',
        action='store',
        dest='covs',
        default='covariates.tsv.gz',
        help='''
        TSV of additional covariates to use. Rows shold be samples. Columns
        should be numeric covariates to include. (default: %(default)s).
        '''
    )

    parser.add_argument(
        '-pf', '--peer_factors',
        action='store',
        dest='pf',
        default='',
        help='''
        TSV containing PEER factor information. Rows should be PEER factors and
        columns should be samples.
        '''
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='merged_covariates.tsv.gz',
        help='Output file.'
    )
    options = parser.parse_args()

    # Get data
    covs = pd.read_csv(
        options.covariates,
        sep='\t',
        header=0,
        index_col=0
    )
    peer_factors = pd.read_csv(
        options.pf,
        sep='\t',
        header=0,
        index_col=0
    ).T
    peer_factors = peer_factors[
        [x for x in peer_factors.columns if x.startswith('factor_')]
    ]

    # Combine
    final_df = pd.concat([covs, peer_factors], axis=1)

    # Do a quick check
    if final_df.shape[0] != covs.shape[0]:
        sys.exit((
            "The number of samples in the final dataframe doesn't match the "
            "number in the original covaraite dataframe."
        ))
    if final_df.shape[0] != peer_factors.shape[0]:
        sys.exit((
            "The number of samples in the final dataframe doesn't match the "
            "number in the original peer factors dataframe."
        ))

    final_df.to_csv(
        options.of,
        index=True,
        header=True,
        sep='\t',
        compression='gzip',
    )
    return


if __name__ == '__main__':
    main()
