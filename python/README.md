# Python Scripts

This directory contains several helper scripts for `RNA-SeQC`.

## Plotting Code

The script `plot.py` can be used to generate a report with summary figures for one or more runs of RNA-SeQC.

`metrics.py` and `nb_encode.py` are both dependencies of `plot.py`.

### Requirements

`plot.py` requires the following python packages:

* numpy
* pandas
* matplotlib
* scipy
* pyBigWig
* bx-python
* agutil
* nbformat
* seaborn
* sklearn

### Usage

```
python plot.py [-h] samples [samples ...] output

positional arguments:
samples     Directory path(s) to RNA-SeQC output. Each directory should
            contain the output files from RNA-SeQC for a single sample
output      Output python notebook

optional arguments:
-h, --help  show this help message and exit
```

Plot.py takes a variable number of arguments. The first and subsequent arguments (excluding the final argument) should be paths to directories where RNA-SeQC output files are located. Each directory should contain output friles from a single run of RNA-SeQC. The final argument should be a filepath to an `.ipynb` file where Plot.py will write its output. The output will contain code cells and pre-rendered output of various summary plots, including  KDE plots of the insert size distributions of each sample and gene expression PCAs.

## Insert size distribution

The script `insert_size_intervals.py` generates a list of intervals in BED format for estimating insert size distributions.

### Requirements

`insert_size_intervals.py` requires the following python packages:

* numpy
* pyBigWig
* qtl

### Usage

```
python3 insert_size_intervals.py [-h] [--min_length MIN_LENGTH]
                                 [--min_mappability MIN_MAPPABILITY]
                                 [--output_dir OUTPUT_DIR]
                                 gtf_path mappability_bigwig prefix

positional arguments:
  gtf_path              Reference annotation in GTF format.
  mappability_bigwig    Mappability track in bigWig format.
  prefix                Prefix for output file names.

optional arguments:
  -h, --help            show this help message and exit
  --min_length MIN_LENGTH
                        Minimum exon/UTR length for computing insert sizes.
                        Default: 1000bp
  --min_mappability MIN_MAPPABILITY
                        Minimum mappability for retained intervals. 
                        Default: 0.95
  --output_dir OUTPUT_DIR
                        Output directory.
```
Mappability tracks can be obtained from [UCSC](https://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability) for GRCh37/hg19 and at [gs://gtex-resources/mappability](https://console.cloud.google.com/storage/browser/gtex-resources/mappability) for GRCh38/hg38.

## Legacy Exon Naming

The script `legacy_exon_remap.py` can be used to translate exon numbers into exon numbers that would have been reported by `RNA-SeQC 1.1.9`. It takes a GTF and Exon GCT and modifies the GCT in place to assign legacy compatible exon names.

### Requirements

`legacy_exon_remap.py` requires the following python packages:

* `https://github.com/francois-a/rnaseq-utils`
* numpy
* pandas
* matplotlib
* scipy
* pyBigWig
* bx-python
* agutil

### Usage

```
python legacy_exon_remap.py [-h] gct gtf

positional arguments:
gct         RNA-SeQC 2 Exon reads gct file
gtf         Reference GTF for the exons

optional arguments:
-h, --help  show this help message and exit
```

Legacy_exon_remap.py takes exactly two arguments. The first should be the filepath to an `exon_reads.gct` file from RNA-SeQC 2. The second should be the same GTF used by RNA-SeQC when the given Exon Reads file was produced. The GCT will be modified in-place to assign RNA-SeQC 1.1.9 compliant exon names.
