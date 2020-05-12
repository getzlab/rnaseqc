#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import glob
import gzip
import os


def combine_gcts(path_dict, verbose=True):
    """Aggregate single-sample GCT files."""

    sample_ids = sorted(path_dict)

    # load first sample and determine dtype
    sample_id = sample_ids[0]
    df = pd.read_csv(path_dict[sample_id], sep='\t', skiprows=3, header=None,
                     index_col=0, names=['Name','Description', sample_id])
    if df[sample_id].dtype == np.float64:
        dtype = np.float32
    elif df[sample_id].dtype == np.int64:
        dtype = np.int32
    else:
        dtype = df[sample_id].dtype.type

    # allocate
    gct_df = pd.DataFrame(0, index=df.index, columns=['Description']+list(sample_ids), dtype=dtype)
    gct_df['Description'] = df['Description']
    gct_df[sample_id] = df[sample_id].astype(dtype)

    for k,sample_id in enumerate(sample_ids[1:], 2):
        if verbose:
            print('\r  * loading GCT {}/{}'.format(k, len(path_dict)), end='', flush=True)
        df = pd.read_csv(path_dict[sample_id], sep='\t', skiprows=3, header=None,
                         usecols=[0,2],  index_col=0, names=['Name', sample_id],
                         dtype={'Name':str, sample_id:dtype})
        gct_df[sample_id] = df[sample_id]
    if verbose:
        print()

    return gct_df


def write_gct(df, gct_file, float_format='%.6g', compresslevel=6):
    """Write pd.DataFrame to GCT format"""

    assert df.index.name=='Name' and df.columns[0]=='Description'

    if gct_file.endswith('.gct.gz'):
        opener = gzip.open(gct_file, 'wt', compresslevel=compresslevel)
    else:
        opener = open(gct_file, 'w')

    with opener as gct:
        gct.write('#1.2\n{0:d}\t{1:d}\n'.format(df.shape[0], df.shape[1]-1))
        df.to_csv(gct, sep='\t', float_format=float_format)


def combine_metrics(path_dict):
    """Aggregate single-sample metrics files."""
    metrics_df = []
    for k,sample_id in enumerate(sorted(path_dict), 1):
        metrics_df.append(pd.read_csv(path_dict[sample_id], sep='\t', index_col=0).astype(np.float32))
    metrics_df = pd.concat(metrics_df, axis=1).T
    metrics_df.index.name = 'sample_id'
    return metrics_df


def combine_insert_sizes(path_dict):
    """Aggregate single-sample insert sizes distributions."""
    insertsize_df = []
    for k,sample_id in enumerate(sorted(path_dict), 1):
        insertsize_df.append(pd.read_csv(path_dict[sample_id], sep='\t', index_col=0, squeeze=True).rename(sample_id))
    insertsize_df = pd.concat(insertsize_df, axis=1).fillna(0).astype(np.int32)
    return insertsize_df



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Aggregate RNA-SeQC outputs')
    parser.add_argument('results_dir', help='Directory containing RNA-SeQC outputs for all samples to be combined.')
    parser.add_argument('prefix', help='Prefix for output files, e.g., <prefix>.gct.gz')
    parser.add_argument('--parquet', action='store_true', help='Write to parquet format instead of GCT')
    parser.add_argument('-o', '--output-dir', default='.', help='Output directory')
    args = parser.parse_args()

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)

    # if os.path.isdir(args.results):
    gene_reads_gcts =  {os.path.basename(i).split('.')[0]:i for i in glob.glob(os.path.join(args.results_dir, '**/*gene_reads.gct*'), recursive=True)}
    gene_fragm_gcts =  {os.path.basename(i).split('.')[0]:i for i in glob.glob(os.path.join(args.results_dir, '**/*gene_fragments.gct*'), recursive=True)}
    gene_tpm_gcts =    {os.path.basename(i).split('.')[0]:i for i in glob.glob(os.path.join(args.results_dir, '**/*gene_tpm.gct*'), recursive=True)}
    exon_reads_gcts =  {os.path.basename(i).split('.')[0]:i for i in glob.glob(os.path.join(args.results_dir, '**/*exon_reads.gct*'), recursive=True)}
    metrics_files =    {os.path.basename(i).split('.')[0]:i for i in glob.glob(os.path.join(args.results_dir, '**/*metrics.tsv*'), recursive=True)}
    insertsize_files = {os.path.basename(i).split('.')[0]:i for i in glob.glob(os.path.join(args.results_dir, '**/*fragmentSizes.txt*'), recursive=True)}
    # coverage files don't get aggregated
    # coverage_files =   {os.path.basename(i).split('.')[0]:i for i in glob.glob(os.path.join(args.results, '*coverage.tsv*'))}

    if len(metrics_files) > 0:
        print('Aggregating metrics')
        metrics_df = combine_metrics(metrics_files)
        metrics_df.to_csv(os.path.join(args.output_dir, '{}.metrics.txt.gz'.format(args.prefix)), sep='\t', float_format='%.6g')

    if len(insertsize_files) > 0:
        print('Aggregating insert size distributions')
        insertsize_df = combine_insert_sizes(insertsize_files)
        insertsize_df.to_csv(os.path.join(args.output_dir, '{}.insert_size_hists.txt.gz'.format(args.prefix)), sep='\t')

    if len(gene_reads_gcts) > 0:
        print('Aggregating read count GCTs')
        gct_df = combine_gcts(gene_reads_gcts)
        if args.parquet:
            gct_df.to_parquet(os.path.join(args.output_dir, '{}.gene_reads.parquet'.format(args.prefix)))
        else:
            write_gct(gct_df, os.path.join(args.output_dir, '{}.gene_reads.gct.gz'.format(args.prefix)))

    if len(gene_fragm_gcts) > 0:
        print('Aggregating fragment count GCTs')
        gct_df = combine_gcts(gene_fragm_gcts)
        if args.parquet:
            gct_df.to_parquet(os.path.join(args.output_dir, '{}.gene_fragments.parquet'.format(args.prefix)))
        else:
            write_gct(gct_df, os.path.join(args.output_dir, '{}.gene_fragments.gct.gz'.format(args.prefix)))

    if len(gene_tpm_gcts) > 0:
        print('Aggregating TPM GCTs')
        gct_df = combine_gcts(gene_tpm_gcts)
        if args.parquet:
            gct_df.to_parquet(os.path.join(args.output_dir, '{}.gene_tpm.parquet'.format(args.prefix)))
        else:
            write_gct(gct_df, os.path.join(args.output_dir, '{}.gene_tpm.gct.gz'.format(args.prefix)))

    if len(exon_reads_gcts) > 0:
        print('Aggregating exon read count GCTs')
        gct_df = combine_gcts(exon_reads_gcts)
        if args.parquet:
            gct_df.to_parquet(os.path.join(args.output_dir, '{}.exon_reads.parquet'.format(args.prefix)))
        else:
            write_gct(gct_df, os.path.join(args.output_dir, '{}.exon_reads.gct.gz'.format(args.prefix)))
