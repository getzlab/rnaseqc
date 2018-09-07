import matplotlib
matplotlib.use("Agg")
import nb_encode as nbe
import argparse
import os
import glob
import sys

# expression = pd.concat(
#     [
#         pd.read_csv(glob.glob(sample+'/*gene_tpm.gct')[0], sep='\t', header=2, index_col=0)[['TPM']].rename({'TPM':sample}, axis='columns')
#         for sample in samples.index
#     ],
#     axis=1
# )
# stats_df = pd.concat([
#         pd.read_csv(path, sep='\t', index_col=0).T.assign(sample_id=[os.path.basename(os.path.dirname(path))]) for path in glob('./*/*.metrics.tsv')
# ]).reset_index().set_index('sample_id')

#
# metrics = pd.concat([
#     pd.read_csv(path+'.metrics.tsv', sep='\t', index_col=0).T.assign(
#         sample_id=sample
#     )
#     for sample, path in samples.items()
# ]).reset_index().set_index('sample_id')
# metrics.columns.name = None


def main(args):
    samples = []
    # TODO: Detect output settings
    # 1) Check RPKM vs TPM in filenames
    # 2) Check for presence of coverage file
    # 3) Check for presence of fragmentSizes file
    for sample in args.samples:
        if not os.path.isdir(sample):
            sys.exit("'%s' is not a valid directory" % sample)
        if not len(glob.glob(sample+"/*.metrics.tsv")):
            sys.exit("'%s' does not contain valid RNA-SeQC output files" % sample)
        samples.append(os.path.join(
            os.path.abspath(sample),
            '.metrics.tsv'.join(
                os.path.basename(glob.glob(sample+"/*.metrics.tsv")[0]).split('.metrics.tsv')[:-1]
            )
        ))
    samples = {os.path.basename(sample):sample for sample in samples}
    print("Detected", len(samples), "samples")

    nb = nbe.Notebook()
    nb.add_markdown_cell(
        '# RNA-SeQC Output',
        '---',
        'Created automatically by RNA-SeQC using the nb_encode api'
    )
    # --- from here, each cell will be run in the script and encoded to the notebook
    import matplotlib.pyplot as plt
    import seaborn as sea
    import pandas as pd
    import numpy as np
    import metrics as met
    nb.add_code_cell([
        'import matplotlib.pyplot as plt',
        'import seaborn as sea',
        'import pandas as pd',
        'import numpy as np',
        'import sys',
        'sys.path.append("%s")' % os.path.dirname(os.path.abspath(__file__)),
        'import metrics as met',
        'import os'
    ])
    # ---
    print("Loading metrics")
    metrics = pd.concat([
        pd.read_csv(path+'.metrics.tsv', sep='\t', index_col=0).T
        for sample, path in samples.items()
    ])
    metrics.columns.name = None
    metrics.index.name = 'sample_id'
    nb.add_code_cell([
        'samples = %s' % repr(samples),
        nbe.trim("""metrics = pd.concat([
            pd.read_csv(path+'.metrics.tsv', sep='\\t', index_col=0).T
            for sample, path in samples.items()
        ])"""),
        'metrics.columns.name = None',
        'metrics.index.name = "sample_id"',
        'metrics.head()'
    ], nbe.encode_dataframe(metrics.head(), nb.exec_count))
    # ---
    print("Generating Fragment Size KDE")
    fig = plt.figure(figsize=(15,10))
    ax = fig.add_subplot(111)
    for sample, path in samples.items():
        path = path+'.fragmentSizes.txt'
        if os.path.isfile(path):
            with open(path) as r:
                line = r.readline().strip()
                if 'Count' in line:
                    # V2
                    fragments = [
                        int(line.strip().split()[0])
                        for line in r
                        for _ in int(line.strip().split()[1])
                    ]
                else:
                    r.seek(0,0)
                    fragments = [int(line.strip()) for line in r]
                if len(np.unique(fragments)) > 1:
                    lim = np.percentile(fragments, 99)
                    sea.kdeplot([x for x in fragments if x <= lim], label=sample, ax=ax)
    ax.set_xlabel("Fragment Length")
    ax.set_ylabel("Density")
    nb.add_code_cell(
        nbe.trim("""
        fig = plt.figure(figsize=(15,10))
        ax = fig.add_subplot(111)
        for sample, path in samples.items():
            path = path+'.fragmentSizes.txt'
            if os.path.isfile(path):
                with open(path) as r:
                    line = r.readline().strip()
                    if 'Count' in line:
                        # V2
                        fragments = [
                            int(line.strip().split()[0])
                            for line in r
                            for _ in int(line.strip().split()[1])
                        ]
                    else:
                        r.seek(0,0)
                        fragments = [int(line.strip()) for line in r]
                    if len(np.unique(fragments)) > 1:
                        lim = np.percentile(fragments, 99)
                        sea.kdeplot([x for x in fragments if x <= lim], label=sample, ax=ax)
        ax.set_xlabel("Fragment Length")
        ax.set_ylabel("Density")
        """),
        ax,
        nbe.encode_figure(fig)
    )
    # ---
    print("Generating QC figures")
    figures = met.plot_qc_figures(
        metrics,
        mapping_rate=0.95,
        million_mapped_reads=40,
        rrna_rate=0.15,
        end1_mismatch_rate=0.005,
        end2_mismatch_rate=0.02,
        intergenic_rate=0.05,
        exonic_rate=0.7,
        alpha=0.5,
        ms=16
    )
    nb.add_code_cell(
        nbe.trim("""
        figures = met.plot_qc_figures(
            metrics,
            mapping_rate=0.95,
            million_mapped_reads=40,
            rrna_rate=0.15,
            end1_mismatch_rate=0.005,
            end2_mismatch_rate=0.02,
            intergenic_rate=0.05,
            exonic_rate=0.7,
            alpha=0.5,
            ms=16
        )
        """),
        *[
            nbe.encode_figure(fig)
            for fig, ax in figures
        ],
        metadata={'scrolled': False}
    )
    # ---
    print("Generating expression PCA")
    expression_df = pd.concat(
        [
            pd.read_csv(path+'.gene_reads.gct', sep='\t', header=2, index_col=0)[['Counts']].rename({'Counts':sample}, axis='columns')
            for sample, path in samples.items()
        ],
        axis=1
    )
    nb.add_code_cell(
        nbe.trim("""
        expression_df = pd.concat(
            [
                pd.read_csv(path+'.gene_reads.gct', sep='\\t', header=2, index_col=0)[['Counts']].rename({'Counts':sample}, axis='columns')
                for sample, path in samples.items()
            ],
            axis=1
        )
        expression_df.head()
        """),
        nbe.encode_dataframe(
            expression_df.head(),
            nb.exec_count
        )
    )
    # ---
    p_df, pca = met.get_pcs(expression_df)
    fig = met.plot_pca(p_df, pca)
    nb.add_code_cell(
        [
            'p_df, pca = met.get_pcs(expression_df)',
            'fig = met.plot_pca(p_df, pca)'
        ],
        nbe.encode_figure(fig),
        metadata={'scrolled': False}
    )
    # ===
    nb.add_code_cell('')
    nb.write(args.output)



if __name__ == '__main__':
    parser = argparse.ArgumentParser('rnaseqc-plot')
    parser.add_argument(
        'samples',
        nargs='+',
        help="Directory path(s) to RNA-SeQC output. Each directory should"
        " contain the output files from RNA-SeQC for a single sample"
    )
    parser.add_argument(
        'output',
        type=argparse.FileType('w'),
        help="Output python notebook"
    )
    main(parser.parse_args())
