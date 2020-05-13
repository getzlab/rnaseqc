# Author: Francois Aguet
import numpy as np
import os
import pyBigWig
import argparse
import qtl.annotation as annotation


def intersect_overlap(intervals):
    """
    intervals: list of tuples or 2-element lists

    breaks intersections into separate intervals
    e.g.: [0,6],[2,8] ->[0,1],[2,6],[7,8]
    """
    intervals = intervals.copy()
    intervals.sort(key=lambda x: (x[0],x[1]))
    intersected = []
    union = list(intervals[0])
    bounds = [intervals[0]]
    for i in intervals[1:]:
        if i[0] <= union[1]:  # overlap w/ previous
            if i[1] > union[1]:  # only extend if larger
                union[1] = i[1]
            bounds.append(i)
        else:
            # process bounds
            if len(bounds)>1:
                p = np.unique([i[0] for i in bounds]+[i[1]+1 for i in bounds])
                intersected.extend([[i,j-1] for i,j in zip(p[:-1],p[1:])])
            else:
                intersected.append(bounds[0])
            # reset
            bounds = [i]
            union = list(i)
    # process last
    if len(bounds)>1:
        p = np.unique([i[0] for i in bounds]+[i[1]+1 for i in bounds])
        intersected.extend([[i,j-1] for i,j in zip(p[:-1], p[1:])])
    else:
        intersected.append(bounds[0])

    return intersected


def parse_intervals(annot, mappability_bw, output_dir, prefix, min_length=1000, min_mappability=0.95):
    """Write intervals to BED format"""

    exclude = set(['retained_intron', 'readthrough_transcript'])

    bw = pyBigWig.open(mappability_bw)
    gintervals = {}
    for c in annot.chr_list:
        exon_coords = []
        for g in annot.chr_genes[c]:
            for t in g.transcripts:
                if (t.type!='retained_intron') and (('tags' not in t.attributes) or len(set(t.attributes['tags']).intersection(exclude))==0):
                    for e in t.exons:
                        exon_coords.append((e.start_pos, e.end_pos))

        v = np.array(intersect_overlap(exon_coords))  # intersect all exons on current chr
        l = v[:,1]-v[:,0]+1
        gintervals[c] = v[l >= min_length, :]
        # filter by mappability
        gintervals[c] = np.array([i for i in gintervals[c] if bw.stats(c, int(i[0])-1, int(i[1]), exact=True)[0] >= min_mappability])
    bw.close()

    # all intervals
    with open(os.path.join(output_dir, '{}_geq{}bp.bed'.format(prefix, min_length)), 'w') as f:
        f.write('#chr\tstart\tend\n')
        for c in annot.chr_list:
            for i in range(gintervals[c].shape[0]):
                f.write('{0:s}\t{1:d}\t{2:d}\n'.format(c, gintervals[c][i][0]-1, gintervals[c][i][1]))  # BED is 0-indexed, [..)

    # single-isoform genes
    gintervals_1iso = {}
    for c in annot.chr_list:
        ec = []
        for g in annot.chr_genes[c]:
            if len(g.transcripts) == 1:
                for e in g.transcripts[0].exons:
                    if e.length >= min_length:
                        ec.append([e.start_pos, e.end_pos])
        ec = list(set([tuple(i) for i in ec]).intersection(set([tuple(i) for i in gintervals[c]])))
        ec.sort(key=lambda x: (x[0],x[1]))
        gintervals_1iso[c] = np.array(ec)

    with open(os.path.join(output_dir, '{}_geq{}bp_1iso.bed'.format(prefix, min_length)), 'w') as f:
        f.write('#chr\tstart\tend\n')
        for c in annot.chr_list:
            for i in range(gintervals_1iso[c].shape[0]):
                f.write('{0:s}\t{1:d}\t{2:d}\n'.format(c, gintervals_1iso[c][i][0]-1, gintervals_1iso[c][i][1]))


if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Parse long exons/UTRs with high mappability for estimating insert size distribution.')
    parser.add_argument('gtf_path', help='Reference annotation in GTF format.')
    parser.add_argument('mappability_bigwig', help='Mappability track in bigWig format.')
    parser.add_argument('prefix', help='Prefix for output file names.')
    parser.add_argument('--min-length', type=np.int32, default=1000, help='Minimum exon/UTR length for computing insert sizes. Default: 1000bp')
    parser.add_argument('--min-mappability', type=np.float64, default=0.95, help='Minimum mappability for retained intervals. Default: 0.95')
    parser.add_argument('--output-dir', default='.', help='Output directory.')
    args = parser.parse_args()

    annot = annotation.Annotation(args.gtf_path, verbose=True)
    parse_intervals(annot, args.mappability_bigwig, args.output_dir, args.prefix,
                    min_length=args.min_length, min_mappability=args.min_mappability)
