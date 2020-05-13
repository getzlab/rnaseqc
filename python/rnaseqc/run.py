#!/usr/bin/env python3
# Author: Francois Aguet
import argparse
import subprocess
from datetime import datetime
import os
import sys

def locate_rnaseqc():
    """
    Search PATH for executable, then try up two directories
    (where compiled executable would exist in the git repo)
    """
    for path in os.environ['PATH'].split(os.pathsep):
        exe = os.path.join(path, 'rnaseqc')
        if test_rnaseqc(exe):
            return exe
    exe = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'rnaseqc')
    if test_rnaseqc(exe):
        return exe
    print("Unable to find rnaseqc executable", file=sys.stderr)
    return 'rnaseqc' # Just try it and see what happens, I guess

def test_rnaseqc(path):
    return os.path.isfile(path) and os.access(path, os.X_OK) and subprocess.run([path, '--version'], stdout=subprocess.PIPE).stdout.startswith(b'RNASeQC 2')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Wrapper for RNA-SeQC 2')
    parser.add_argument('genes_gtf', type=str, help='Gene annotation GTF')
    parser.add_argument('bam_file', type=str, help='BAM file')
    parser.add_argument('prefix', type=str, default='Reads', help='Prefix for output files; usually sample_id')
    parser.add_argument('-o', '--output_dir', default=os.getcwd(), help='Output directory')
    parser.add_argument('-q', '--mapping-quality', default=None, type=int, help="Lower bound on read quality for reads used in coverage metrics")
    parser.add_argument('-m', '--mismatch-threshold', default=None, type=int, help="Maximum allowed mismatches in a read while still used for coverage metrics")
    parser.add_argument('-c', '--coverage', action='store_true', help="Include raw coverage metrics in a separate output table. By default, only summary statistics are included in metrics")
    parser.add_argument('--stranded', default=None, choices=['rf', 'fr'], help='Strandedness for stranded libraries')
    parser.add_argument('--bed', default=None, help='BED file with intervals for estimating insert size distribution')
    args = parser.parse_args()

    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running RNA-SeQC', flush=True)

    cmd = '{} {} {} {}'.format(locate_rnaseqc(), args.genes_gtf, args.bam_file, args.output_dir) \
        + ' -s '+args.prefix \
        + ' -vv'
    if args.stranded is not None:
        cmd += ' --stranded '+args.stranded
    if args.bed is not None:
        cmd += ' --bed '+args.bed
    if args.mapping_quality is not None:
        cmd += ' --mapping-quality {}'.format(args.mapping_quality)
    if args.mismatch_threshold is not None:
        cmd += ' --base-mismatch {}'.format(args.mismatch_threshold)
    if args.coverage:
        cmd += ' --coverage'
    print('  * command: "{}"'.format(cmd), flush=True)
    subprocess.check_call(cmd, shell=True)

    # gzip GCTs
    subprocess.check_call('gzip {0}.exon_reads.gct {0}.gene_tpm.gct {0}.gene_reads.gct {0}.gene_fragments.gct'.format(args.prefix), shell=True)
    if args.coverage:
        subprocess.check_call('gzip {}.coverage.tsv'.format(args.prefix), shell=True)

    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finished RNA-SeQC', flush=True)
