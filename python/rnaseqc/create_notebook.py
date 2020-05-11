from . import nb_encode as nbe
import argparse
import subprocess


def main(args):

    nb = nbe.Notebook()

    nb.add_markdown_cell(
        '# RNA-SeQC metrics report',
    )

    nb.add_code_cell([
        'import pandas as pd',
        'import qtl.io',
        'import rnaseqc.report'
    ])

    nb.add_code_cell([
        "# load inputs",
        "metrics_df = pd.read_csv('{}', sep='\\t', index_col=0)".format(args.metrics),
        "tpm_df = qtl.io.read_gct('{}')".format(args.tpm),
        "cohort_s = pd.read_csv('{}', sep='\\t', index_col=0, header=None, squeeze=True)".format(args.cohort),
        "insertsize_df = pd.read_csv('{}', sep='\\t', index_col=0)".format(args.insert_size),
    ])

    nb.add_code_cell([
        "thresholds = {'Exonic Rate': 0.7}",
        'rnaseqc.report.plot_qc_figures(metrics_df, cohort_s=cohort_s, cohort_colors=None, date_s=None,',
        '                               show_legend=True, ms=12, alpha=1, highlight_ids=None,',
        '                               thresholds=thresholds, insertsize_df=insertsize_df, tpm_df=tpm_df)',
    ])

    nb.add_code_cell('')
    nb.write(args.output)



if __name__ == '__main__':
    parser = argparse.ArgumentParser('rnaseqc-plot')
    parser.add_argument('metrics', help='Aggregated metrics')
    parser.add_argument('output', type=argparse.FileType('w'),
                        help="Output python notebook")
    parser.add_argument('-t', '--tpm', default=None, help='Aggregated TPM')
    parser.add_argument('-i', '--insert-size', default=None,
                        help='Aggregated insert size distributions')
    parser.add_argument('-c', '--cohort', default=None,
                        help='TSV file mapping sample IDs to cohort/batch IDs')
    parser.add_argument('-d', '--date', default=None,
                        help='TSV file mapping sample IDs to dates')
    args = parser.parse_args()

    # generate notebook
    main(args)

    # execute notebook
    subprocess.check_call('jupyter nbconvert --execute --inplace {}'.format(args.output.name), shell=True)
