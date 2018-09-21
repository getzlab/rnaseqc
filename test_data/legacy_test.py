import argparse
import subprocess
import pandas as pd
import numpy as np


def main(args):
    df = pd.read_csv(args.modern, index_col=0, header=2, sep='\t').join(
        pd.read_csv(args.legacy, index_col=0, header=2, sep='\t')[['RNA-SeQC']]
    )
    if args.approx:
        assert len(df[np.abs(df['Counts'] - df['RNA-SeQC']) > .01]) == 0
    else:
        assert len(df[df['Counts'] != df['RNA-SeQC']]) == 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser('legacy-test')
    parser.add_argument(
        'modern',
        help='Output file from modern tool'
    )
    parser.add_argument(
        'legacy',
        help='Output file from legacy (java) tool'
    )
    parser.add_argument(
        '-a', '--approx',
        action='store_true',
        help="Check that values are within .01 of eachother. Good for exon counts"
        " which usually differ by java rounding precision"
    )
    args = parser.parse_args()
    main(args)
