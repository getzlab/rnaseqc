import argparse
import subprocess
import pandas as pd
import numpy as np


def main(args):
    if args.mode == 'tables':
        df = pd.read_csv(args.input1, index_col=0, header=2, sep='\t').join(
            pd.read_csv(args.input2, index_col=0, header=2, sep='\t'),
            how='outer',
            rsuffix='_'
        )
    elif args.mode == 'metrics' or args.mode == 'fragments':
        df = pd.read_csv(args.input1, sep='\t', index_col=0).join(
            pd.read_csv(args.input2, sep='\t', index_col=0),
            how='outer',
            rsuffix='_',

        )
    assert len(df[np.abs(df[args.columns[0]] - df[args.columns[1]]) > args.tolerance]) == 0, df[np.abs(df[args.columns[0]] - df[args.columns[1]]) > args.tolerance].head()
    if args.mode == 'fragments':
        assert len(set(pd.read_csv(args.input1, sep='\t', index_col=0).index) ^ set(pd.read_csv(args.input2, sep='\t', index_col=0).index)) == 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser('legacy-test')
    parser.add_argument(
        'input1',
        help='First input file'
    )
    parser.add_argument(
        'input2',
        help='Second input file'
    )
    parser.add_argument(
        '-t', '--tolerance',
        nargs='?',
        type=float,
        help="Tolerance for differing values. If not provided, this defaults to 0,"
        " for exact comparison. If provided without argument, this defaults to .01,"
        " which is usually good for checking modern vs legacy counts (which vary"
        " slightly within Java's default precision). You can also provide a floating"
        " point number to manually specify tolerance",
        default=0.000001,
        const=0.01
    )
    parser.add_argument(
        '-m', '--mode',
        choices=[
            'metrics',
            'tables',
            'fragments'
        ],
        default='metrics',
        help="What type of input file is being compared. Default: metrics"
    )
    parser.add_argument(
        '-c', '--columns',
        nargs=2,
        help="Column names to load for 'tables'",
        metavar=['COLUMN-A', 'COLUMN-B'],
        default=('Counts', 'RNA-SeQC')
    )
    args = parser.parse_args()
    main(args)
