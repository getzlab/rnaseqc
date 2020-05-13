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
        'report': get_filepath('report'),
        'run': get_filepath('run'),
    }

    run_parser = subparsers.add_parser(
        'run',
        help='A light wrapper with some convenience functions to run RNA-SeQC',
        description='A light wrapper with some convenience functions to run RNA-SeQC',
        add_help=False
    )

    aggregate_parser = subparsers.add_parser(
        'aggregate',
        help='Aggregate RNA-SeQC outputs from multiple samples',
        description='Aggregate RNA-SeQC outputs from multiple samples',
        add_help=False
    )

    notebook_parser = subparsers.add_parser(
        'notebook',
        help='Generate a notebook with figures comparing outputs from multiple samples',
        description='Generate a notebook with figures comparing outputs from multiple samples',
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
        help='Generate a BED file with long (>1000bp), high-mappability intervals for estimating insert sizes',
        description='Generate a BED file with long (>1000bp), high-mappability intervals for estimating insert sizes',
        add_help=False
    )

    legacy_parser = subparsers.add_parser(
        'legacy-exons',
        help='Renames exons in exon_reads.gct file from RNA-SeQC 2 to use naming convention from RNA-SeQC 1.1.x',
        description='Renames exons in exon_reads.gct file from RNA-SeQC 2 to use naming convention from RNA-SeQC 1.1.x',
        add_help=False
    )

    args, remainder = parser.parse_known_args()
    if args.command in commands:
        os.execvp(sys.executable, [sys.executable, '-m', commands[args.command]] + remainder)
    else:
        parser.print_usage()
        sys.exit('A valid subcommand must be provided.')
