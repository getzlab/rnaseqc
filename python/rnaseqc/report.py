import pandas as pd
import numpy as np
import argparse
import os
import sys
import matplotlib.pyplot as plt
import qtl.io

sys.path.insert(1, os.path.dirname(__file__))
from .plot import *


def plot_qc_figures(metrics_df, cohort_s=None, cohort_colors=None, date_s=None,
                    prefix=None, output_dir=None, dpi=300, show_legend=True,
                    ms=12, alpha=1, show_xticklabels=False, highlight_ids=None,
                    thresholds=None, insertsize_df=None, tpm_df=None):
    """
    metrics_df: output from RNA-SeQC
    cohort_s: mapping of sample ID to cohort/cluster/etc.
    """
    if cohort_s is None:
        cohort_s = pd.Series('All samples', index=metrics_df.index)

    if output_dir is not None:
        assert prefix is not None

    cohorts = np.unique(cohort_s)
    if cohort_colors is None:
        cohort_colors = get_cohort_colors(cohorts)

    metrics_args = {
        'cohort_s': cohort_s,
        'cohort_colors': cohort_colors,
        'show_xticklabels': show_xticklabels,
        'ms': ms,
        'alpha': alpha,
        'highlight_ids': highlight_ids
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
        "Median 3' bias": None,
        'Median Exon CV': None,
        'Average Fragment Length': None,
    }

    # plot cohort legend
    if show_legend:
        ax = qtl.plot.setup_figure(0.5,0.5)
        for c in cohorts:
            ax.scatter([], [], s=48, marker='s', c=cohort_colors[c].reshape(1,-1), label=c)
        ax.legend(loc='upper left', handlelength=1, ncol=5, bbox_to_anchor=(1,1))
        plt.axis('off')

    # distributions for selected/key metrics
    for k,metric in enumerate(metrics_list, 1):
        if metric in metrics_df:
            ylim = [0,1]
            metrics(metrics_df[metric], ylim=ylim_dict[metric],
                         threshold=threshold_dict.get(metric, None),
                         threshold_dir=threshold_dir_dict.get(metric, None),
                         **metrics_args)
            if output_dir is not None:
                plt.savefig(os.path.join(output_dir, '{}.{}.pdf'.format(prefix, metric.lower().replace("3'",'3prime').replace(' ','_'))), dpi=dpi)

    # genes detected vs bias and duplication rate
    detection_bias(metrics_df, bias_metric="Median 3' bias", c='Duplicate Rate of Mapped')
    if output_dir is not None:
        plt.savefig(os.path.join(output_dir, '{}.genes_detected_vs_median_3prime_bias.pdf'.format(prefix)), dpi=dpi)

    # mismatch rates
    mismatch_rates(metrics_df, cohort_s=cohort_s, cohort_colors=cohort_colors)
    if output_dir is not None:
        plt.savefig(os.path.join(output_dir, '{}.end_mismatch_rates.pdf'.format(prefix)), dpi=dpi)

    mapping_sense(metrics_df, cohort_s, cohort_colors=cohort_colors)
    if output_dir is not None:
        plt.savefig(os.path.join(output_dir, '{}.mapping_sense.pdf'.format(prefix)), dpi=dpi)

    # insert size distributions (if supplied)
    if insertsize_df is not None:
        insert_sizes(insertsize_df, cohort_s, cohort_colors=cohort_colors, sort_order='mean')
        if output_dir is not None:
            plt.savefig(os.path.join(output_dir, '{}.insert_sizes.pdf'.format(prefix)), dpi=dpi)

    if tpm_df is not None:
        cdf_df = calculate_expression_cdfs(tpm_df)
        if tpm_df.shape[1] < 50:
            mode = 'lines'
        else:
            mode = 'ci'
        cumulative_expression(cdf_df, cohort_s=cohort_s, cohort_colors=cohort_colors, mode=mode)
        if output_dir is not None:
            plt.savefig(os.path.join(output_dir, '{}.cumulative_expression.pdf'.format(prefix)), dpi=dpi)


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
        cohort_s = pd.read_csv(args.cohort, sep='\t', index_col=0, header=None, squeeze=True)
        assert metrics_df.index.isin(cohort_s.index).all()
    else:
        cohort_s = None

    if args.cohort is not None:
        insertsize_df = pd.read_csv(args.insert_size, sep='\t', index_col=0)
    else:
        insertsize_df = None

    return metrics_df, tpm_df, cohort_s, insertsize_df


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate QC report from RNA-SeQC metrics table.')
    parser.add_argument('metrics', help='Aggregated QC metrics from RNA-SeQC.')
    parser.add_argument('prefix', help='Name for output files.')
    parser.add_argument('--tpm', default=None, help='Aggregated TPM matrix from RNA-SeQC.')
    parser.add_argument('--insert-size', default=None, help='Aggregated insert sizes from RNA-SeQC.')
    parser.add_argument('--cohort', default=None, help='Cohort or batch annotation. TSV file mapping sample IDs to annotation.')
    parser.add_argument('--output-dir', default='.', help='If specified, figures are saved here.')
    parser.add_argument('--dpi', type=int, default=300, help='Figure resolution.')
    args = parser.parse_args()

    metrics_df, tpm_df, cohort_s, insertsize_df = load_inputs(args)

    plot_qc_figures(metrics_df, cohort_s=cohort_s, cohort_colors=None, date_s=None,
                    prefix=args.prefix, output_dir=args.output_dir, dpi=args.dpi, show_legend=True,
                    ms=12, alpha=1, show_xticklabels=False, highlight_ids=None,
                    thresholds=None, insertsize_df=insertsize_df, tpm_df=tpm_df)
