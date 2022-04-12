#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import pandas as pd
import numpy as np
import plotnine as plt9


def var_mean_expr_plot(
    df,
    x_axis='mean_expr',
    y_axis='cis',
    fit_formula=None
):
    plt = (
        plt9.ggplot(df, plt9.aes(x=x_axis, y=y_axis))
        + plt9.geom_point(alpha=0.35)
        + plt9.labs(
            x='Mean expression',
            y='Variance explained by {} effects'.format(y_axis)
        )
        + plt9.theme_bw()
    )

    if fit_formula is not None:
        plt = plt + plt9.stat_smooth(method='lm', formula=fit_formula)

    return plt


def mean_expr_quantiles(
    df,
    x_axis='mean_expr',
    y_axis='cis',
    n_bins=10,
    y_max=100,
    y_min=0
):
    df['binned_mean_expr'] = pd.qcut(
        df[x_axis],
        q=n_bins,
        labels=[i+1 for i in range(n_bins)]
    )

    summary_df = df.groupby('binned_mean_expr').agg(
        mean_y=pd.NamedAgg(
            column=y_axis,
            aggfunc=lambda x: np.mean(x)
        ),
        sd_y=pd.NamedAgg(
            column=y_axis,
            aggfunc=lambda x: np.std(x)
        )
    ).reset_index()
    summary_df['max_y'] = summary_df.apply(
        lambda x: min(x['mean_y'] + x['sd_y'], y_max),
        axis=1
    )
    summary_df['min_y'] = summary_df.apply(
        lambda x: max(x['mean_y'] - x['sd_y'], y_min),
        axis=1
    )

    plt = (
        plt9.ggplot(summary_df, plt9.aes(
            x='binned_mean_expr',
            y='mean_y'
        ))
        + plt9.geom_line()
        + plt9.geom_point(alpha=0.75)
        + plt9.geom_errorbar(
            plt9.aes(
                ymax='max_y',
                ymin='min_y'
            ),
            width=.2,
            position="dodge"
        )
        + plt9.theme_bw()
        + plt9.theme(legend_position=None)
        + plt9.labs(
            x='Mean expression bin',
            y='Variance explained by {} effects'.format(y_axis)
        )
        + plt9.ylim(y_min, y_max)
    )
    return(plt)


def plot_effect(
    df,
    effect,
    base_dir
):
    plt = var_mean_expr_plot(
        df,
        'mean_expr',
        effect,
        fit_formula='y ~ x'
    )
    plt.save(
        filename='{}/{}_mean_expr.png'.format(base_dir, effect),
        dpi=320
    )

    # plt = var_mean_expr_plot(
    #     df,
    #     'mean_expr',
    #     effect,
    #     log2_x=True,
    #     log10_y=False,
    #     fit_formula='y ~ x + x**2'
    # )
    # plt.save(
    #     filename='{}/quad_fit__{}_log2_mean_expr.png'.format(
    #       base_dir, effect
    #     ),
    #     dpi=320
    # )

    plt = mean_expr_quantiles(
        df,
        x_axis='mean_expr',
        y_axis=effect,
        n_bins=10,
        y_max=100,
        y_min=0
    )
    plt.save(
        filename='{}/{}_binned_mean_expr.png'.format(base_dir, effect),
        dpi=320
    )


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="Plot variance decomposition results."
    )

    parser.add_argument(
        '-vd', '--vardec_file',
        action='store',
        dest='vd',
        required=True,
        help='''
        TSV file containing variance decomposition information. Should have
        the following columns: noise, mean_expr, and those specified  in
        `effect_labels`.
        '''
    )

    parser.add_argument(
        '-el', '--effect_labels',
        action='store',
        dest='el',
        required=True,
        help='''
        Comma-separated list of effects considered.
        '''
    )

    parser.add_argument(
        '-opd', '--output_plot_dir',
        action='store',
        dest='opd',
        default='plots',
        help='Output directory for plots'
    )
    options = parser.parse_args()

    vardec = pd.read_csv(options.vd, header=0, sep='\t')
    effects = options.el.split(',')

    # First normalize variances
    # See: https://github.com/limix/limix/blob/ebcd261a95ae7e0a3ecb7f17b7d379cc4140b328/limix/vardec/_vardec.py#L214
    var_subset = vardec.loc[:, effects + ['noise']].copy()
    vardec.loc[:, effects + ['noise']] = (
        var_subset.divide(var_subset.sum(axis=1), axis=0) * 100
    )

    for effect in effects:
        _ = plot_effect(vardec, effect, options.opd)


if __name__ == '__main__':
    main()
