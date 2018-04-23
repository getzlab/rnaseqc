# RNA-SeQC
Version 2.0.0

---

**NOTE**: It is recommended that the provided GTF be collapsed in such a way that there are no overlapping transcripts **on the same strand**, otherwise reads may be double-counted in areas where the read aligns to multiple exons on the same strand. This is **not** a transcript-quantification method. Readcounts and coverage are made towards exons and genes only if *all* aligned segments of a read fully align to exons of a gene, but keep in mind that coverage may be counted towards multiple transcripts (and its exons) if these criteria are met. Beyond this, no attempt will be made to disambiguate which transcript a read belongs to.

### Command Line Usage:

`rnaseqc [OPTIONS] gtf bam output`

###### OPTIONS:
      -h, --help                        Display this message and quit

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

      --read-length=[LENGTH]            Set the maximum accepted length. Reads
                                        longer than this threshold are
                                        discarded. Default: 1000000 [bp]

      --fragment-samples=[SAMPLES]      Set the number of samples to take when
                                        computing fragment sizes. Requires the
                                        --bed argument. Default: 1000000

      --low-quality=[QUALITY]           Set the lower bound on read quality.
                                        Reads below this number are counted as
                                        low quality BUT ARE STILL USED IN
                                        COUNTS. See --mapping-quality to discard
                                        reads based on quality. Default: 255

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

      --chimeric-tag=[TAG]              Reads maked with the specified tag will
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
{sample}.metrics.tsv contains all metrics recorded by RNASeQC:
* Mapping Rate: Proportion of all reads in the Bam which are primary alignments, mapped, and have not failed vendor QC (based on flags).
* Unique Rate of Mapped: Proportion of all mapped reads (see above) which are not marked as duplicates (based on flags).
* Duplicate Rate of Mapped: Proportion of all mapped reads (see "Mapping Rate", above) which are marked as duplicates (based on flags). This should be 1 - (Unique Rate of Mapped).
* Base Mismatch: Proportion of bases in mapped reads (see "Mapping Rate", above) which are mismatched (based on the value of the NM tag).
* End 1/2 Mapping Rate: Proportion of mapped reads (see "Mapping Rate", above) which are paired and are the first/second (respectively) mate in the pair.
* End 1/2 Mismatch Rate: Proportion of bases of all first/second mate reads (see above) which are mismatched (based on the value of the NM tag).
* Expression Profiling Efficiency: Proportion of all reads in the bam which are considered "Exonic" (see "Exonic Rate", below).
* Exonic Rate: Proportion of all mapped reads (see "Mapping Rate", above) where all aligned segments of the reads fully aligned to exons of the same gene(s).
* Intronic Rate: Proportion of all mapped reads (see "Mapping Rate", above) where the read at least partially aligned to a gene, but did not align to any exons of the gene(s).
* Intergenic Rate: Proportion of all mapped reads (see "Mapping Rate", above) where the read never aligned to any part of a gene.
* Intragenic Rate: Proportion of all mapped reads (see "Mapping Rate", above) where the read at least partially aligned to a gene. This should equal to (Exonic Rate) + (Intronic Rate), but does not complement the Intergenic Rate (see "Disqualification Rate" and "Discard Rate", below).
* Disqualification Rate: Proportion of all mapped reads (see "Mapping Rate", above) where at least one aligned segment of a read fully aligned to at least one exon, and at least one other aligned segment was intronic, intergenic, or aligned to an exon of a different gene. These reads are disqualified because it is ambiguous which gene the read belongs to.
* Discard Rate: Proportion of all mapped reads (see "Mapping Rate", above) which were discarded due to preliminary filters (the read was longer than the maximum **--read-length**, was considered "chimeric" (see "Chimeric Reads", below) and the **--exclude-chimeric** flag was set, had more mismatches than the set **--base-mismatch** threshold, the **--unpaired** flag was not set and the read was not in a proper pair, the read had a mapping quality lower than the set **--mapping-quality** threshold, or the read posessed any **--tag**s set to be filtered) and never checked against the GTF. The (Exonic Rate), (Intronic Rate), (Intergenic Rate), (Disqualification Rate) and (Discard Rate) should sum to exactly 1
* rRNA Rate: Proportion of all mapped reads (see "Mapping Rate", above) which passed preliminary filters (see above) and where at least one aligned segment of the read at least partially intersected a feature in the GTF where the _transcript\_type_ was set to either "Mt_rRNA" or "rRNA".
* End 1/2 Sense Rate: Proportion of all mapped reads (see "Mapping Rate", above) which passed prelimiary filters (see "Discard Rate", above) and where the read was the first/second mate (respectively) and was aligned to the sense strand (based on if the read aligned to a forward/reverse strand feature, was the first/second mate, and if the read was marked as being on the reverse strand or not).
* Avg. Splits per Read: Average number of indels or introns in each mapped read (see "Mapping Rate", above) which separate the aligned segments of the read.
* Several raw couns are recorded, which are rather self-explanetory or whose meaning can be inferred from the above rates. Select raw counts are explained here:
    * Chimeric Reads: Raw count of mapped reads (see "Mapping Rate", above) which either posessed the **--chimeric-tag** (which defaults to "mC") or where the read was paried, the mate was marked as mapped (by flag), and the mates were either on different contigs or were separated by more than the **--chimeric-distance**. The presense of the chimeric tag overrides categorizing chimera by other conditions. If at least one read posesses the chimeric tag, then the other conditions are ignored. In **--legacy** mode, the chimeric tag is ignored and, in addition to the criteria described above, reads on a contig which appears 127th or below in the Bam's sequence header are also counted as chimeric.
    * Genes Detected: Raw count of genes present in the GTF where at least 5 reads had all of their aligned segments fully align to exons of this gene.
    * Read Length: The maximum read length encountered, excluding any reads which were longer than the maximum **--read-length**
    * Estimated Library Complexity: Library complexity estimated as in PicardTools
* 3' Bias statistics (mean, median, standard deviation, median absolute deviation, and the 25th and 75th percentile of these rates): An estimation of bias towards the 3' end of genes based on Intergenic Reads (see "Intergenic Rate", above) which align (at least partially) to windows at the 3' and 5' ends of each gene. These windows are **--window-size** bases wide, and are offset (inwards) **--offset** bases into the gene, away from the true 3' and 5' ends. Bias calculations are only performed on genes which are at least **--gene-length** bases long.
* Fragment Length statistics (mean, median, standard deviation, and median absolute deviation of fragment lengths): An estimation of fragment lengths (also called insert size) for each pair. Fragment length calculations are done on reads which pass preliminary filters (see "Discard Rate", above) and are computed against the **--bed** file, if is given. Iff all of the aligned segments of a **both** mates in a pair intersect the same interval (and only that interval) of the Bed file, record that pair's insert size. Only the first **--fragment-samples** pairs meeting these criteria will be recorded. The Bed file should contain sufficiently long, non-overlapping, exon and UTR intervals. Fragment size statistics are not reported without a Bed file.
* Transcript Coverage statistics (overall medians of the mean, standard deviation, and coefficient of variance of per-base coverage of each transcript): Aligned segments of reads are used to record per-base coverage of exons from Exonic Reads (see "Exonic Rate", above). As per the note at the top, it is recommended that your GTF contain only one transcript per gene and no overlapping transcripts on the same strand, as it may cause a read to count coverage towards multiple transcripts.
* Delta CV: The normalized median of absolute difference of log2 CV between adjacent exons. A log2 CV is computed for each exon which had coverage, and the absolute difference between all adjacent exons are recorded (never comparing exons from different transcripts).

### Legacy mode differences

The **--legacy** flag enables compatibility with RNASeQC 1.1.9. This ensures that exon and gene readcounts match exactly the counts which would have been produced by running that version. This also adds an extra condition to classify reads as chimeric (see "Chimeric Reads", above). Any metrics which existed in 1.1.9 will also match within Java's floating point precision.

### Code Documentation:

###### Included source files:
* RNASeQC.cpp : Main entry point for the application
* Expression.cpp : Function definitions for running Exon and Gene expression metrics (the bulk of the program)
* Expression.h : Function declarations for Expression.h
* Metrics.cpp : Function definitions for the `Collector`, `Metrics`,  `CoverageEntry`, `BaseCoverage`, and `BiasCounter` classes
* Metrics.h :  `Collector`, `Metrics`,  `CoverageEntry`, `BaseCoverage`, and `BiasCounter` class declarations
* GTF.cpp : Function definitions for parsing GTF files and for `Feature` class functions
* GTF.h : `Feature` class and related function declarations
* BED.cpp : `extractBED()` function definition
* BED.h : Declares `extractBED()`

#### Documentation

##### RNASeQC
* `int main(int, char**)` : Parses command line arguments then runs RNASeQC.  First it reads and parses the GTF into lists of features (one list per contig).  This accounts for the bulk of the memory overhead.  If a BED file was provided, it will also parse it and store it in similar lists.  Then it iterates over each read in the BAM file, and counts various attributes of the reads.  Reads which pass filtering (must be a primary alignment, be mapped, have an overall read span less than the maximum read length, be in a proper pair, have a mapping quality >= the quality threshold, have no more than the maximum number of mismatches, and cannot fail vendor QC) are run through exon coverage metrics (see Expression.cpp).  If a BED file was provided, each read which passes the same filtering steps is also run through fragment size metrics.  After all reads have been processed, some statistics are caluclated and then a report is generated

* `bool compGenes(const string&, const string&)`: Compares genes using provided gene_ids as strings. Returns true iff the first gene has a lower TPM than the second.

* `template <typename T> double computeMedian(unsigned long, T&&, unsigned int = 0u)`: Utility function to compute the median of a sorted container by providing its size and an iterator to its start position (which satisfies `Forward Iterator` requirements).

* `void add_range(vector<unsigned long>&, coord, unsigned int)`: Utility function to increment readcounts in a vector representing the per-base coverage of an exon.

* `tuple<double, double, double> computeCoverage(ofstream&, string&, string&, const double, map<string, vector<CoverageEntry> >&, vector<string>&)`: Function which computes the coverage of a single transcript and returns it's mean, standard deviation, and coefficient of variation as a tuple. If the ofstream is open, it will also write these statistics to a file.

##### Expression
* `unsigned int extractBlocks(BamTools::BamAlignment&, std::vector<Feature>&, unsigned short);` : Function to extract the aligned segments of a read into separate objects for querying the GTF

* `std::list<Feature>* intersectBlock(Feature&, std::list<Feature>&);` : Function to query one read segment against the list of GTF features.  **Note: The returned list is dynamically allocated and must be freed**

* `void trimFeatures(BamTools::BamAlignment&, std::list<Feature> &);` : Function to trim out features from the GTF which will not be queried again. This helps to speed up the program as it advances down each contig

* `unsigned int fragmentSizeMetrics(<args...>)` : Queries each segment of a read against the intervals in the BED file.  If all segments of a read as well as all segments of the read's mate map to the same interval, the pair is deemed 'valid' and its insert size is recorded in the list of fragment sizes

* `void exonAlignmentMetrics(<args...>)` : Queries each segment of a read against the features in the GTF. To be counted towards a gene, all segments of a read must fully map to exons of the gene. Exon and gene coverage is discarded for reads which do not meet that criteria (having bases overhang off an exon/gene, aligning to introns or intergenic regions, etc). Reads are counted as Intronic if none of their segments aligned to any exons (even if the coverage is discarded), and Intergenic if they additionally never intersected any genes. Because of this criteria, a read is counted as disqualified if not all of its segments fully mapped to exons of the same gene(s).

* `void legacyExonAlignmentMetrics(<args...>)` : Modified version of `exonAlignmentMetrics()` which produces counts that match the output of
RNASeQC 1.1.9.

##### Metrics
* `Metrics` (class) : Used for keeping track of various counters and statistics regarding the reads that have been processed so far

* `Collector` (class) : Used to collect exon coverage per gene, per read. Only exons belonging to genes that are not disqualified (all segments of a read map to an exon of the same gene) will be kept.

* `CoverageEntry` (struct) : Used to record a single, contiguous coverage interval of a read on an exon

* `BaseCoverage` (class) : Used to keep track of coverage as it's recorded during execution. Provisional coverage is added to each exon that a read aligns to, and is only comitted if the read is considered "Exonic" (see "Exonic Rate", above) and only for genes where every aligned segment to the read aligns to an exon of that gene. When an exon is trimmed out of the search window (see `void trimFeatures`, in "Expression", above) the coverage of that exon is dumped to a temporary file for later. This helps greatly reduce the memory footprint (usually keeping it below 1Gb) during execution at the cost of a slight performance loss due to the IO. At the end of the program, coverage is read back in from the temporary file and used to compute per-base coverage metrics for each transcript.

* `BiasCounter` (class) : Used to record 3'/5' coverage on each gene and compute a bias ratio based on these coverage levels.

* `std::ofstream& operator<<(std::ofstream&, const Metrics&)` : Dumps all the raw counts from the `Metrics` counter to the output stream

##### GTF
* `coord` (typedef) : `long long`

* `Feature` (struct) : Intended to represent GTF features. However, throughout the program it's used to represent any type of genomic interval

* `bool operator==(const Feature &a, const Feature &b)` : returns true if the two features have the same start, end, chromosome, strand, type, feature id, and transcript id

* `bool compIntervalStart(const Feature&, const Feature&)` : returns true if the first feature starts before the second

* `bool compIntervalEnd(const Feature&, const Feature&)` : returns true if the first feature ends before the second

* `bool intersectPoint(const Feature&, const coord)` : returns true if the coordinate falls within the feature's interval

* `bool intersectInterval(const Feature&, const Feature&)` : returns true if the two intervals have any overlap

* `int partialIntersect(const Feature&, const Feature&)` : returns the exact size of the overlap between two features

* `std::map<std::string, unsigned short> chromosomes` (external variable) : Externally available map of chromosomes to their internal ID.  IDs are arbitrarily assigned (and are used in place of strings for performance reasons)

* `std::map<std::string, std::string> geneNames` (external variable) : Externally available map of gene IDs to their canonical names.  Used for generating gct output files

* `std::map<std::string, coord> geneLengths` (external variable) : Externally available map of gene IDs to their reference length.  Used for generating TPM and RPKM values

* `std::map<std::string, coord> transcriptCodingLengths` (external variable) : Externally available map of gene IDs to their reference length.  Used for generating TPM and RPKM values

* `std::vector<std::string> geneList` (external variable) : List of gene ids

* `std::vector<std::string> exonList` (external variable) : List of exon ids

* `unsigned short chromosomeMap(std::string)` : returns the internal ID for the provided chromosome name

* `std::ifstream& operator>>(std::ifstream&, Feature&)` : reads from the ifstream (assumes GTF format).  Stores the next feature encountered in the provided reference

* `std::map<std::string,std::string>& parseAttributes(std::string&, std::map<std::string,std::string>&)` : Parses the attribute string of a GTF feature into a map

##### BED
* `std::ifstream& extractBED(std::ifstream&, Feature&)` : reads from the ifstream (assumes BED format). Stores the next feature encountered in the provided reference
