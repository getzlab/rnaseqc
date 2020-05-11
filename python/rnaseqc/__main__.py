import argparse
import sys
import os
import subprocess

def get_filepath(name):
    # return os.path.join(os.path.dirname(os.path.abspath(__file__)), name)
    return 'rnaseqc.{}'.format(name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser('rnaseqc')

    subparsers = parser.add_subparsers(dest='command')

    commands = {
        'aggregate': get_filepath('aggregate'),
        'notebook': get_filepath('create_notebook'),
        'insert-size': get_filepath('insert_size_intervals'),
        'legacy-exons': get_filepath('legacy_exon_remap'),
        'report': get_filepath('report')
    }

    aggregate_parser = subparsers.add_parser(
        'aggregate',
        help='Aggregate RNA-SeQC outputs from multiple samples',
        description='Aggregate RNA-SeQC outputs from multiple samples',
        add_help=False
    )

    notebook_parser = subparsers.add_parser(
        'notebook',
        help='Generate a notebook with figures to compare outputs from multiple samples',
        description='Generate a notebook with figures to compare outputs from multiple samples',
        add_help=False
    )

    report_parser = subparsers.add_parser(
        'report',
        help='Generate PDF figures from aggregated RNA-SeQC results',
        description='Generate PDF figures from aggregated RNA-SeQC results',
        add_help=False
    )

    intervals_parser = subparsers.add_parser(
        'insert-size',
        help='Generate a BED file containing "good" intervals for RNA-SeQC\'s insert-size estimation metrics',
        description='Generate a BED file containing "good" intervals for RNA-SeQC\'s insert-size estimation metrics',
        add_help=False
    )

    legacy_parser = subparsers.add_parser(
        'legacy-exons',
        help='Renames exons in exon_reads.gct file from RNA-SeQC 2.X.X to use exon names as reported by RNA-SeQC 1.1.9',
        description='Renames exons in exon_reads.gct file from RNA-SeQC 2.X.X to use exon names as reported by RNA-SeQC 1.1.9',
        add_help=False
    )

    args, remainder = parser.parse_known_args()
    if args.command in commands:
        os.execvp(sys.executable, [sys.executable, '-m', commands[args.command]] + remainder)
    else:
        parser.print_usage()
        sys.exit("You must provide a valid subcommand")
