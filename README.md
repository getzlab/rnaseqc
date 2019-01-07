# RNA-SeQC
[![Build Status](https://travis-ci.com/broadinstitute/rnaseqc.svg?token=y8NpD4Ye6EkYyigZUZDt&branch=master)](https://travis-ci.com/broadinstitute/rnaseqc)

Version 2.1.0

## Installing

The latest stable build of RNA-SeQC is available on the [GitHub Releases](https://github.com/broadinstitute/rnaseqc/releases) page, and contains static binaries for Linux and OSX.

RNA-SeQC is also available as a docker image: `gcr.io/broad-cga-aarong-gtex/rnaseqc:latest` which is automatically updated with any code change.
Older versions of the docker image are tagged using the full commit SHA of any commit which introduced a code change.

To checkout the source of RNA-SeQC run `git clone --recursive https://github.com/broadinstitute/rnaseqc.git`.
If you do not use the `--recursive` flag, you'll need to run `git submodule update --init --recursive` or you will be missing [SeqLib](https://github.com/walaj/SeqLib).

#### Unit Tests

If you checkout the RNA-SeQC source and wish to run unit tests, you'll need to have [Git LFS](https://git-lfs.github.com/) installed.
If you already have LFS installed when you clone the repository, it will automatically
pull the test resources. If you install LFS after the fact, just run `git lfs pull` to
manually download the test data. The test resources use **~850 MB** of space.

You can run the unit tests with `make test`

## Usage

**NOTE**: This tool requires that the provided GTF be collapsed in such a way that there are no overlapping transcripts **on the same strand** and that each gene have a single transcript whose id matches the parent gene id. This is **not** a transcript-quantification method. Readcounts and coverage are made towards exons and genes only if *all* aligned segments of a read fully align to exons of a gene, but keep in mind that coverage may be counted towards multiple transcripts (and its exons) if these criteria are met. Beyond this, no attempt will be made to disambiguate which transcript a read belongs to.
You can collapse an existing GTF using the [GTEx collapse annotation script](https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model)

### Command Line Usage:

`rnaseqc [OPTIONS] gtf bam output`

Example: `./rnaseqc test_data/downsampled.gtf test_data/downsampled.bam --bed test_data/downsampled.bed --coverage .`

###### OPTIONS:
      -h, --help                        Display this message and quit

      --version                         Display the version and quit

      gtf                               The input GTF file containing features
                                        to check the bam against

      bam                               The input SAM/BAM file containing reads
                                        to process

      output                            Output directory

      -s[sample], --sample=[sample]     The name of the current sample. Default:
                                        The bam's filename

      --bed=[BEDFILE]                   Optional input BED file containing
                                        non-overlapping exons used for fragment
                                        size calculations

      --chimeric-distance=[DISTANCE]    Set the maximum accepted distance
                                        between read mates. Mates beyond this
                                        distance will be counted as chimeric
                                        pairs. Default: 2000000 [bp]

      --fragment-samples=[SAMPLES]      Set the number of samples to take when
                                        computing fragment sizes. Requires the
                                        --bed argument. Default: 1000000

      -q[QUALITY],
      --mapping-quality=[QUALITY]       Set the lower bound on read quality for
                                        exon coverage counting. Reads below this
                                        number are excluded from coverage
                                        metrics. Default: 255

      --base-mismatch=[MISMATCHES]      Set the maximum number of allowed
                                        mismatches between a read and the
                                        reference sequence. Reads with more than
                                        this number of mismatches are excluded
                                        from coverage metrics. Default: 6

      --split-distance=[DISTANCE]       Set the minimum distance between aligned
                                        blocks of a read for the read to be counted
                                        as split. Default: 100 [bp]

      --offset=[OFFSET]                 Set the offset into the gene for the 3'
                                        and 5' windows in bias calculation. A
                                        positive value shifts the 3' and 5'
                                        windows towards eachother, while a
                                        negative value shifts them apart.
                                        Default: 150 [bp]

      --window-size=[SIZE]              Set the size of the 3' and 5' windows in
                                        bias calculation. Default: 100 [bp]

      --gene-length=[LENGTH]            Set the minimum size of a gene for bias
                                        calculation. Genes below this size are
                                        ignored in the calculation. Default:
                                        600 [bp]

      --legacy                          Use legacy read counting rules. Read
                                        counts match output of RNA-SeQC 1.1.9

      --stranded=[stranded]             Use strand-specific metrics. Only
                                        features on the same strand of a read
                                        will be considered. Allowed values are
                                        'RF', 'rf', 'FR', and 'fr'

      -v, --verbose                     Give some feedback about what's going
                                        on. Supply this argument twice for
                                        progress updates while parsing the bam

      -t[TAG...], --tag=[TAG...]        Filter out reads with the specified tag.

      --chimeric-tag=[TAG]              Reads marked with the specified tag will
                                        be labeled as Chimeric. Defaults to 'mC'
                                        for compatibility with STAR

      --exclude-chimeric                Exclude chimeric reads from the read
                                        counts

      -u, --unpaired                    Treat all reads as unpaired, ignoring
                                        filters which require properly paired
                                        reads

      --rpkm                            Output gene RPKM values instead of TPMs

      --coverage                        If this flag is provided, coverage
                                        statistics for each transcript will be
                                        written to a table. Otherwise, only
                                        summary coverage statistics are
                                        generated and added to the metrics table

      --coverage-mask=[SIZE]            Sets how many bases at both ends of a
                                        transcript are masked out when computing
                                        per-base exon coverage. Default: 500bp

      -d[threshold],
      --detection-threshold=[threshold] Number of counts on a gene to consider
                                        the gene 'detected'. Additionally, genes
                                        below this limit are excluded from 3'
                                        bias computation. Default: 5 reads

      "--" can be used to terminate flag options and force all following
      arguments to be treated as positional options

### Output files:
The following output files are generated in the output directory you provide:
* {sample}.metrics.tsv : A tab-delimited list of (Statistic, Value) pairs of all statistics and metrics recorded.
* {sample}.exon_reads.gct : A tab-delimited GCT file with (Exon ID, Gene Name, coverage) tuples for all exons which had at least part of one read mapped.
* {sample}.gene_reads.gct : A tab-delimited GCT file with (Gene ID, Gene Name, coverage) tuples for all genes which had at least one read map to at least one of its exons
* {sample}.gene_tpm.gct : A tab-delimited GCT file with (Gene ID, Gene Name, TPM) tuples for all genes reported in the gene_reads.gct file. Note: this file is renamed to .gene_rpkm.gct if the **--rpkm** flag is present.
* {sample}.fragmentSizes.txt : A list of fragment sizes recorded, if a BED file was provided
* {sample}.coverage.tsv : A tab-delimited list of (Gene ID, Transcript ID, Mean Coverage, Coverage Std, Coverage CV) tuples for all transcripts encountered in the GTF.

#### Metrics reported:

See [Metrics.md](Metrics.md) for a description of all metrics reported in the `metrics.tsv`, `coverage.tsv`, and `fragmentSizes.txt` files.

### Legacy mode differences

The **--legacy** flag enables compatibility with RNASeQC 1.1.9. This ensures that exon and gene readcounts match exactly the counts which would have been produced by running that version. This also adds an extra condition to classify reads as chimeric (see "Chimeric Reads", above). Any metrics which existed in 1.1.9 will also match within Java's floating point precision.
