#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import pandas as pd
import numpy as np


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="Prepare moltraits file for VarDec."
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
        '-ss', '--sample_start',
        action='store',
        dest='ss',
        type=int,
        required=True,
        help='''
        Column index specifying start of sample data. 0-indexed.
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
        default='moltraits.tsv',
        help='Output file.'
    )
    options = parser.parse_args()

    # Grab data
    moltraits = pd.read_csv(options.mt, sep='\t', header=0)
    samples = [
        str(x) for x in # in case the samples are numbers. Pandas will read headers as string
        pd.read_csv(options.samples, sep='\n', header=None)[0].to_list()
    ]

    # First check to make sure all samples are found in moltraits
    missing_samp = np.setdiff1d(samples, moltraits.columns.to_list()).tolist()
    if len(missing_samp) > 0:
        raise ValueError(
            'Some samples specified in {} are not in {}: {}'.format(
                options.samples,
                options.mt,
                ', '.join(missing_samp)
            )
        )

    #
    cols_to_retain = moltraits.columns[0:options.ss].tolist()
    cols_to_retain.extend(samples)
    filtered_moltraits = moltraits[cols_to_retain].copy()

    filtered_moltraits.to_csv(
        options.of,
        header=True,
        index=False,
        sep='\t',
        compression='gzip'
    )


if __name__ == '__main__':
    main()
