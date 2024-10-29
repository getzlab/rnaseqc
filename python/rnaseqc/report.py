import pandas as pd
import numpy as np
import argparse
import os
import sys
import matplotlib.pyplot as plt
import qtl.io

sys.path.insert(1, os.path.dirname(__file__))
from .plot import *


def plot_qc_figures(metrics_df, cohort_s=None, cohort_order=None, cohort_colors=None, date_s=None,
                    insertsize_df=None, gc_content_df=None, tpm_df=None,
                    thresholds=None, lims=None, outlier_method='threshold',
                    show_legend=True, legend_cols=5, lw=4, lh=1, ms=12, alpha=1, show_xticklabels=False,
                    highlight_ids=None, prefix=None, output_dir=None, dpi=300):
    """
    metrics_df: output from RNA-SeQC
    cohort_s: mapping of sample ID to cohort/cluster/etc.
    """
    if cohort_s is None:
        cohort_s = pd.Series('All samples', index=metrics_df.index)
    else:
        assert metrics_df.index.isin(cohort_s.index).all() and cohort_s.loc[metrics_df.index].notnull().all()

    if date_s is not None:
        assert metrics_df.index.isin(date_s.index).all() and date_s.loc[metrics_df.index].notnull().all()

    if output_dir is not None:
        assert prefix is not None

    cohorts = np.unique(cohort_s.loc[metrics_df.index])
    if cohort_colors is None:
        cohort_colors = get_cohort_colors(cohorts)

    metrics_args = {
        'cohort_s': cohort_s,
        'cohort_order': cohort_order,
        'cohort_colors': cohort_colors,
        'date_s': date_s,
        'show_xticklabels': show_xticklabels,
        'ms': ms,
        'alpha': alpha,
        'highlight_ids': highlight_ids,
        'aw': 6,
        'ah': 2,
    }

    metrics_list = [
        'Mapped Reads',  # Unique Mapping, Vendor QC Passed Reads that were mapped
        'Mapping Rate',
        'Duplicate Rate of Mapped',
        'Exonic Rate',
        'Intronic Rate',
        'Intergenic Rate',
        'Chimeric Alignment Rate',
        'rRNA Rate',
        # 'Mapped Unique Reads',  # Duplicate Rate of Mapped is more representative
        "Median 3' bias",
        'Median Exon CV',
        'Fragment GC Content Mean',
        'Average Fragment Length',
    ]

    threshold_dir_dict = {
        'Mapped Reads': 'lt',
        'Mapping Rate': 'lt',
        'Duplicate Rate of Mapped': 'gt',
        'Exonic Rate': 'lt',
        'Intronic Rate': 'gt',
        'Intergenic Rate': 'gt',
        'Chimeric Alignment Rate': 'gt',
        'rRNA Rate': 'gt',
        "Median 3' bias": 'gt',
        'Median Exon CV': 'gt',
        'Average Fragment Length': 'lt',
    }

    threshold_dict = {
        'Mapped Reads': 50e6,
        'Mapping Rate': 0.9,
        # 'Duplicate Rate of Mapped': 0.6,
        'Exonic Rate': 0.7,
        'Intronic Rate': 0.05,
        'Intergenic Rate': 0.1,
        'Chimeric Alignment Rate': 0.01,
        'rRNA Rate': 0.1,
    }
    if thresholds is not None:
        threshold_dict.update(thresholds)

    ylim_dict = {
        'Mapped Reads': None,
        'Mapping Rate': [0,1],
        'Duplicate Rate of Mapped': [0, 1],
        'Exonic Rate': [0, 1],
        'Intronic Rate': [0, 1],
        'Intergenic Rate': [0, 1],
        'Chimeric Alignment Rate': [0, 0.1],
        'rRNA Rate': [0, 1],
        "Median 3' bias": [0, 1],
        'Median Exon CV': None,
        'Fragment GC Content Mean': [0, 1],
        'Average Fragment Length': None,
    }

    if lims is not None:
        ylim_dict.update(lims)

    if cohort_order is None:
        cohort_order = cohorts

    # plot cohort legend
    if show_legend:
        ax = qtl.plot.setup_figure(lw, lh, xspace=[0,0], yspace=[0,0])
        for c in cohort_order:
            ax.scatter(np.nan, np.nan, s=48, marker='s', color=cohort_colors[c], label=c)
        ax.scatter(np.nan, np.nan, fc='w', ec='k', lw=1, s=30, label='Outliers')
        ax.legend(loc='center left', handlelength=1, ncol=legend_cols)
        plt.axis('off')
        if output_dir is not None:
            plt.savefig(os.path.join(output_dir, f'{prefix}.legend.pdf'), dpi=dpi)

    # distributions for selected/key metrics
    for k,metric in enumerate(metrics_list, 1):
        if metric in metrics_df and not (metrics_df[metric] == 0).all():
            if metric == 'Duplicate Rate of Mapped' and 'Duplicate Rate of Mapped, excluding Globins' in metrics_df:
                metric_s = metrics_df['Duplicate Rate of Mapped, excluding Globins'].rename('Duplicate Rate of Mapped')
            else:
                metric_s = metrics_df[metric]
            metrics(metric_s, ylim=ylim_dict[metric],
                         threshold=threshold_dict.get(metric, None),
                         threshold_dir=threshold_dir_dict.get(metric, None),
                         outlier_method=outlier_method,
                         **metrics_args)
            if output_dir is not None:
                plt.savefig(os.path.join(output_dir, '{}.{}.pdf'.format(prefix, metric.lower().replace("3'",'3prime').replace(' ','_'))), dpi=dpi)

    # genes detected vs bias and duplication rate
    if "Median 3' bias" in metrics_df:
        c = 'Duplicate Rate of Mapped, excluding Globins' if 'Duplicate Rate of Mapped, excluding Globins' in metrics_df else 'Duplicate Rate of Mapped'
        detection_bias(metrics_df, bias_metric="Median 3' bias", c=c)
        if output_dir is not None:
            plt.savefig(os.path.join(output_dir, f'{prefix}.genes_detected_vs_median_3prime_bias.pdf'), dpi=dpi)

    # mismatch rates
    if not metrics_df['End 1 Mismatch Rate'].isnull().all():
        mismatch_rates(metrics_df, cohort_s=cohort_s, cohort_order=cohort_order, cohort_colors=cohort_colors,
                       end1_threshold=threshold_dict.get('End 1 mismatch rate', None),
                       end2_threshold=threshold_dict.get('End 2 mismatch rate', None))
        if output_dir is not None:
            plt.savefig(os.path.join(output_dir, f'{prefix}.end_mismatch_rates.pdf'), dpi=dpi)

    mapping_sense(metrics_df, cohort_s=cohort_s, cohort_order=cohort_order,
                  cohort_colors=cohort_colors, date_s=date_s, width=1)
    if output_dir is not None:
        plt.savefig(os.path.join(output_dir, f'{prefix}.mapping_sense.pdf'), dpi=dpi)

    # insert size distributions (if supplied)
    if insertsize_df is not None:
        insert_sizes(insertsize_df, cohort_s=cohort_s, cohort_order=cohort_order,
                     cohort_colors=cohort_colors, sort_order='cohort')
        if output_dir is not None:
            plt.savefig(os.path.join(output_dir, f'{prefix}.insert_sizes.pdf'), dpi=dpi)

    if gc_content_df is not None:
        gc_content(gc_content_df, cohort_s=cohort_s, cohort_colors=cohort_colors,
                   cohort_order=cohort_order, sort_order='cohort')
        if output_dir is not None:
            plt.savefig(os.path.join(output_dir, f'{prefix}.gc_content.pdf'), dpi=dpi)

    if tpm_df is not None:
        cdf_df = calculate_expression_cdfs(tpm_df)
        if tpm_df.shape[1] < 50:
            mode = 'lines'
        else:
            mode = 'ci'
        cumulative_expression(cdf_df, cohort_s=cohort_s, cohort_colors=cohort_colors, mode=mode)
        if output_dir is not None:
            plt.savefig(os.path.join(output_dir, f'{prefix}.cumulative_expression.pdf'), dpi=dpi)


def load_inputs(args):

    if args.metrics.endswith('.parquet'):
        metrics_df = pd.read_parquet(args.metrics)
    else:
        metrics_df = pd.read_csv(args.metrics, sep='\t', index_col=0)

    if args.tpm is not None:
        tpm_df = qtl.io.read_gct(args.tpm, load_description=False)
    else:
        tpm_df = None

    if args.cohort is not None:
        cohort_s = pd.read_csv(args.cohort, sep='\t', index_col=0, header=None).squeeze('columns')
        assert metrics_df.index.isin(cohort_s.index).all()
    else:
        cohort_s = None

    if args.date is not None:
        date_s = pd.read_csv(args.date, sep='\t', index_col=0, header=None).squeeze('columns')
        assert metrics_df.index.isin(date_s.index).all()
    else:
        date_s = None

    if args.insert_size is not None:
        insertsize_df = pd.read_csv(args.insert_size, sep='\t', index_col=0)
    else:
        insertsize_df = None

    return metrics_df, tpm_df, cohort_s, date_s, insertsize_df


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate QC report from RNA-SeQC metrics table.')
    parser.add_argument('metrics', help='Aggregated QC metrics from RNA-SeQC.')
    parser.add_argument('prefix', help='Name for output files.')
    parser.add_argument('--tpm', default=None, help='Aggregated TPM matrix from RNA-SeQC.')
    parser.add_argument('--insert-size', default=None, help='Aggregated insert sizes from RNA-SeQC.')
    parser.add_argument('--cohort', default=None, help='Cohort or batch annotation. TSV file mapping sample IDs to annotation.')
    parser.add_argument('--date', default=None, help='Date annotation. TSV file mapping sample IDs to dates.')
    parser.add_argument('--output-dir', default='.', help='If specified, figures are saved here.')
    parser.add_argument('--dpi', type=int, default=300, help='Figure resolution.')
    args = parser.parse_args()

    metrics_df, tpm_df, cohort_s, date_s, insertsize_df = load_inputs(args)

    plot_qc_figures(metrics_df, cohort_s=cohort_s, cohort_colors=None, date_s=date_s,
                    prefix=args.prefix, output_dir=args.output_dir, dpi=args.dpi, show_legend=True,
                    ms=12, alpha=1, show_xticklabels=False, highlight_ids=None,
                    thresholds=None, insertsize_df=insertsize_df, tpm_df=tpm_df)
