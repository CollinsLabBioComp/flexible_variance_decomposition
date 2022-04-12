#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import pandas as pd
import os.path


def mean_compute_columns(df, columns):
    # Iterate through each column, get NaN's and replace with mean
    for col in columns:
        df.loc[pd.isna(df[col]), col] = df[col].mean()

    return df


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="Prepare covariates for VarDec."
    )

    parser.add_argument(
        '-cov_f', '--covariate_file',
        action='store',
        dest='covariate_file',
        default='',
        help='''
        TSV file containing numeric covariates. Rows must be samples and columns
        must be covariates. First column should be sample IDs.
        '''
    )

    parser.add_argument(
        '-covs', '--covariates',
        action='store',
        dest='covs',
        default=None,
        help='''
        Comma-separated list of covariates to retain in downstream models.
        '''
    )

    parser.add_argument(
        '-covs_ref', '--covariate_references',
        action='store',
        dest='covs_ref',
        default=None,
        help='''
        Comma-separated list of covariate levels. For each covariate, it should
        have the format: covariate_1::reference,covariate_2::references.
        '''
    )

    parser.add_argument(
        '-covs_mi', '--covariate_mean_impute',
        action='store',
        dest='covs_mi',
        default=None,
        help='''
        Comma-separated list of covariates to mean impute.
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
        default='covariates.tsv.gz',
        help='Output file.'
    )
    options = parser.parse_args()

    # Grab data
    samples = [
        str(x) for x in # in case the samples are numbers. Pandas will read headers as string
        pd.read_csv(options.samples, sep='\n', header=None)[0].to_list()
    ]

    if os.path.isfile(options.covariate_file) and options.covs is not None:
        covs = pd.read_csv(
            options.covariate_file,
            sep='\t',
            header=0,
            index_col=0
        )
        covs = covs.loc[samples, ]  # grab samples

        # Now limit covariates to specified covs
        covs_to_retain = options.covs.split(',')
        covs = covs[covs_to_retain]

        # Get dictonary for reference values
        ref_dict = {}
        if options.covs_ref is not None:
            for ref_cov in options.covs_ref.split(','):
                tmp_val = ref_cov.split('::')
                ref_dict[tmp_val[0]] = tmp_val[1]


        for col in covs.columns:
            if col in ref_dict.keys():
                cats = [ref_dict[col]]
                _ = [cats.append(x) for x in covs[col] if x not in cats]
                covs[col] = [cats.index(x) for x in covs[col]]

                # Print mapping
                pstr = ['{}: {}'.format(cats.index(x), x) for x in cats]
                print('Casting column {} to continuous. Index mapping:'.format(
                    col
                ))
                print('\n'.join(pstr))
            # Cast any that are characters -- both PEER and LIMIX require cont.
            elif not pd.api.types.is_numeric_dtype(covs[col]):
                unique_vals = covs[col].unique().tolist()
                covs[col] = [unique_vals.index(x) for x in covs[col]]

                # Print mapping
                pstr = ['{}: {}'.format(
                    unique_vals.index(x), x
                ) for x in unique_vals]
                print('Casting column {} to continuous. Index mapping:'.format(
                    col
                ))
                print('\n'.join(pstr))

        # Mean compute columns if need be
        if options.covs_mi is not None:
            mean_comp_cols = options.covs_mi.split(',')
            covs = mean_compute_columns(covs, mean_comp_cols)
    else:
        # Create blank dataframe from samples
        covs = pd.DataFrame(index=samples)

    # Add 1's for offset
    # Write as 1.0 for automatic conversion to double in R
    covs['offset'] = [1.0] * covs.shape[0]

    covs.to_csv(
        options.of,
        index=True,
        header=True,
        sep='\t',
        compression='gzip',
    )


if __name__ == '__main__':
    main()
