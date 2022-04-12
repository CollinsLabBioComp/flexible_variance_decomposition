#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import pandas as pd


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="Merge list of dataframes."
    )

    parser.add_argument(
        '-dfs', '--dataframes',
        action='store',
        dest='dfs',
        required=True,
        help='''
        Comma-separated list of dataframes to merge. Dataframes should be
        TSVs. Files should contain headers.
        '''
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='merged_df.tsv.gz',
        help='Output file with merged dataframes.'
    )
    options = parser.parse_args()

    # Get data
    files = options.dfs.split(',')

    final_df = pd.read_csv(files[0], sep='\t', header=0)
    for file in files[1:]:
        df = pd.read_csv(file, sep='\t', header=0)
        final_df = final_df.append(df, ignore_index=True)

    final_df.to_csv(
        options.of,
        index=False,
        header=True,
        sep='\t',
        compression='gzip'
    )
    return


if __name__ == '__main__':
    main()
