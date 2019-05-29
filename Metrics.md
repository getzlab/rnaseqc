# RNA-SeQC Output Metrics

This file provides a description for each of the output metrics in the `metrics.tsv` file. A description of other output files can be found at the bottom

## Output Metrics
* Mapping Rate: The proportion of all reads in the Bam which were Mapped, and not Secondary Alignments or Platform/Vendor QC Failing reads ("Mapped Reads").
In legacy mode, this is the proportion of all reads which were Mapped out
of all reads which were not Secondary Alignments or Platform/Vendor QC Failing reads.
* Unique Rate of Mapped: This is the proportion of reads which **were not** marked as PCR/Optical Duplicates out of all "Mapped Reads" (as defined above; excludes Secondary and Vendor QC Failed reads).
* Duplicate Rate of Mapped: This is the proportion of all reads which **were** marked as PCR/Optical Duplicates out of all "Mapped Reads" (as defined above; excludes Secondary and Vendor QC Failed reads). This is complementary to the "Unique Rate of Mapped".
* Duplicate Rate of Mapped, excluding Globins: This is similar to the "Duplicate Rate of Mapped" except that it only includes reads which **did not** align to _HBA1_, _HBA2_, _HBB_, or _HBD_.
* Base Mismatch: The total number of mismatched bases (as determined by the "NM" tag) of all "Mapped Reads" (as defined above) divided by the total aligned length of all "Mapped Reads".
* End 1 & 2 Mapping Rate: The proportion of Paired reads which were marked as First or Second in the pair, respectively, out of all "Mapped Reads" (above).
* End 1 & 2 Mismatch Rate: The proportion of mismatched bases (as determined by the "NM" tag) belonging to First or Second mates, divided by the total aligned length of all "Mapped" (above) First or Second mates, respectively.
* Expression Profiling Efficiency: The proportion of "Exonic Reads" (see "Exonic Rate", below) out of all reads which were not Secondary Alignments or
Platform/Vendor QC Failing reads.
* High Quality Rate: The proportion of **properly paired** reads with less than 6 mismatched bases and a perfect mapping quality out of all "Mapped Reads" (above).
* Exonic Rate: The proportion of "Mapped Reads" (above) for which all aligned segments unambiguously aligned to exons of the same gene.
* Intronic Rate: The proportion of "Mapped Reads" (above) for which all aligned segments unambiguously aligned to the same gene, but none of which _intersected_ any exons of the gene.
* Intergenic Rate: The proportion of "Mapped Reads" (above) for which none of the aligned segments _intersected_ any genes.
* Intragenic Rate: The sum of "Exonic" and "Intronic" rates (see "Exonic Rate" and "Intronic Rate" above)
* Ambiguous Alignment Rate: The proportion of "Mapped Reads" (above) where the aligned segments unambiguously aligned to exons of more than one gene.
* High Quality Exonic, Intronic, Intergenic, Intragenic, and Ambiguous Alignment Rates: The proportion of "Exonic Reads", "Intronic Reads", "Intragenic Reads", "Intergenic Reads", and "Ambiguous Reads" (see rates above) out of "High Quality Reads" only (as defined in "High Quality Rate", above)
* Discard Rate: The proportion of "Mapped Reads" (above) which discarded and not checked against the reference annotation. In most cases this should be 0, however, this will include reads which were discarded by additional command line flags (such as `--exclude-chimeric` or `--tag`) or extra legacy mode filters. "Exonic Rate", "Intronic Rate", "Intergenic Rate", "Ambiguous Alignment Rate" and "Discard Rate" will sum to 1.
* rRNA Rate: The proportion of "Mapped Reads" (above) which at least partially intersected with an annotated rRNA gene. This is **not** complementary to any other rates.
* End 1 & 2 Sense Rate: The proportion of First or Second Mate reads which intersected with a Sense Strand feature out of all First or Second
Mate reads which intersected with any features, respectively.
* Avg. Splits per Read: The average number of gaps or deletions present in "Mapped Reads" (above). This is generally not an important metric, but may indicate aligner errors if the value is too high.
* Raw Counts: All raw counts used in any metrics are reported here in the file
* Read Length: The longest aligned length observed in any read
* Genes Detected: The number of genes which had at least 5 unambiguous reads. The detection threshold can be changed with `--detection-threshold`
* Estimated Library Complexity: An estimation of the number of unique cDNA fragments present in the library. This computation follows the same formula as Picard EstimateLibraryComplexity
* 3' Bias statistics (Mean, Median, Std Deviation, Median Absolute Deviation, 25th percentile, 75th percentile): These aggregate statistics are based on the total coverage in 100 bp windows on both the 3' and 5' ends of a gene. The windows are both offset 150 bases into the gene. This computation is only performed on genes at least 600bp long and with at least 5 unambiguous reads. These thresholds can be changed with `--offset`, `--window-size`, `--gene-length`, and `--detection-threshold`. A gene with even coverage in both it's 3' and 5' windows would have a bias of 0.5; bias near 1 or 0 may indicate degredation
* Fragment Length Statistics (Mean, Meadian, Std Deviation, and Median Absolute Deviation): These aggregate statistics are based on the insert sizes observed in "High Quality" (above) read pairs. These metrics are only present if a Bed file was provided with the `--bed` option. Only the first 1,000,000 "High Quality" pairs, where both mates map to the same Bed interval are used.
* Median of Transcript Coverage statistics (Mean, Std Deviation, Coefficient of Variation): These statistics are the median of a given aggregate statistic of transcript coverage (for example, the median of mean transcript coverage). Transcript coverage is computed by dropping the first and last 500bp of each gene and measuring the **"High Quality"** (above) coverage over the remainder of the gene.
* Median Exon CV: The median coefficient of variation of exon coverage. Exon coverage is computed by dropping the first and last 500bp of each gene and measuring the **"High Quality"** (above) coverage over the remainder of the exons. This is considered a good metric for sample quality. A lower value indicates more consistent coverage over exons.
* Exon CV MAD: The Median Absolute Deviation over all Exon CVs

**Note**: When running in `--unpaired` mode, single-ended bams will report `nan` for all End 1 and End 2 metrics

### Fragment Sizes File

This file contains the raw counts of the observed insert sizes of the sample. Fragment sizes are only measured if a Bed file is provided with the `--bed` option. This file is stored as a histogram, with the first column recording a given observed size, and the second column recording the number of occurances of that particular size.

### Coverage File

This file contains coverage data for all genes. Coverage computations are always performed, but this file of per-gene coverage data is not produced unless
the `--coverage` flag is provided. The first column contains the gene ID as given by the input annotation. The next three columns contain the mean, standard deviation, and coefficient of variation of coverage for each gene, respectively. The first and last 500bp of each gene are dropped and not considered when computing coverage. A value of 0 or `nan` may indicate that the gene's coding length was less than 1kb or that the gene had 0 coverage
over it's exons.

## Migrating between old and new columns

For users of the legacy tool, several metrics have been renamed, removed, or changed.
Below is a table of previous metrics and how to access them using the new metrics names:

Old Metric | New Metric | Notes
-|-|-
Base Mismatch Rate | Base Mismatch |
Duplication Rate of Mapped | Duplicate Rate of Mapped |
End 1/2 % Sense | End 1/2 Sense Rate |
Estimated Library Size | Esitmated Library Complexity |
Failed Vendor QC Check | Failed Vendor QC |
Fragment Length Mean | Average Fragment Length | The fragment length metrics have changed significantly
Fragment Length StdDev | Fragment Length Std
Intragenic Rate | Intragenic Rate | Some reads previously classified as `Intragenic` are now classified as `Ambiguous Alignments`. The equivalent of the old `Intragenic Rate` can be computed by summing `Intragenic Rate` + `Ambigous Alignment Rate`
Mapped | Mapped Reads |
Mapped Unique | Mapped Unique Reads |
Total Purity Filtered Reads Sequenced | Unique Mapping, Vendor QC Passed Reads | This counts reads without the Secondary or QC Fail flags set. For a true count of total alignments use `Total Reads`
