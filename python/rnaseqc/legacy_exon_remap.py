import argparse
from agutil import status_bar
import subprocess
import csv
import shutil
from qtl.annotation import Annotation
import tempfile

def run(args):
    print("Parsing GTF")
    gtf = Annotation(args.gtf.name)
    print("Parsing GCT")
    numRows = int(subprocess.check_output("wc -l %s" % args.gct.name, shell=True).decode().strip().split()[0]) - 3
    header = ''.join([next(args.gct), next(args.gct)])
    reader = csv.DictReader(args.gct, delimiter='\t')
    w = tempfile.NamedTemporaryFile('w')
    w.write(header)
    writer = csv.DictWriter(w, reader.fieldnames, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    current = None
    features = []
    with status_bar(numRows) as bar:
        for line in reader:
            bar.update(bar.current + 1)
            gene = '_'.join(line['Name'].split('_')[:-1])
            if gene != current:
                if current is not None:
                    ref = gtf.get_gene(current)
                    try:
                        if len(ref):
                            ref = ref[0]
                    except:
                        pass
                    exons = {exon.id:exon for transcript in ref.transcripts for exon in transcript.exons}
                    raw_size = len(exons)
                    for exon in [exon for exon in exons]:
                        try:
                            if exon.isdigit() and int(exon) <= raw_size:
                                exons[current+'_'+exon] = exons[exon]
                        except:
                            pass
                    features.sort(
                        key=lambda feat:(
                            1 if exons[feat['Name']].length == 1 else 0,
                            exons[feat['Name']].start_pos,
                            exons[feat['Name']].end_pos
                        )
                    )
                    for i in range(len(features)):
                        parts = features[i]['Name'].split('_')
                        prefix = '_'.join(parts[:-1])
                        suffix = parts[-1]
                        if exons[features[i]['Name']].length == 1:
                            features[i][reader.fieldnames[-1]] = 0
                        suffix = str(i)
                        features[i]['Name'] = prefix+'_'+suffix
                    writer.writerows(features)
                current = gene
                features = []
            features.append({k:v for k,v in line.items()})
    if len(features):
        ref = gtf.get_gene(current)
        try:
            if len(ref):
                ref = ref[0]
        except:
            pass
        exons = {exon.id:exon for transcript in ref.transcripts for exon in transcript.exons}
        raw_size = len(exons)
        for exon in [exon for exon in exons]:
            try:
                if exon.isdigit() and int(exon) <= raw_size:
                    exons[current+'_'+exon] = exons[exon]
            except:
                pass
        features.sort(
            key=lambda feat:(
                1 if exons[feat['Name']].length == 1 else 0,
                exons[feat['Name']].start_pos,
                exons[feat['Name']].end_pos
            )
        )
        for i in range(len(features)):
            prefix, suffix = features[i]['Name'].split('_')
            if exons[features[i]['Name']].length == 1:
                features[i]['Counts'] = 0
            suffix = str(i)
            features[i]['Name'] = prefix+'_'+suffix
            writer.writerows(features)
    print("Cleaning up")
    w.flush()
    args.gct.close()
    shutil.copyfile(
        args.gct.name,
        args.gct.name+'.bak'
    )
    shutil.copyfile(
        w.name,
        args.gct.name
    )


def main():
    parser = argparse.ArgumentParser('flipper')

    parser.add_argument(
        'gct',
        type=argparse.FileType('r'),
        help="RNA-SeQC 2 Exon reads gct file"
    )

    parser.add_argument(
        'gtf',
        type=argparse.FileType('r'),
        help="Reference GTF for the exons"
    )

    args = parser.parse_args()
    run(args)

if __name__ == '__main__':
    main()
