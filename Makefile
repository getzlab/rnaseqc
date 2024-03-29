#Set inclusion paths here (if boost, bamtools, or args are installed outside your path)
INCLUDE_DIRS=-ISeqLib -ISeqLib/htslib/
#Set library paths here (if boost or bamtools are installed outside your path)
LIBRARY_PATHS=
#Set to 0 if you encounter linker errors regarding strings from the bamtools library
ABI=1
#Provide full paths here to .a archives for libraries which should be statically linked
STATIC_LIBS=
#List of remaining libraries that will be dynamically linked
LIBS= -lboost_filesystem -lboost_regex -lboost_system -lz -llzma -lbz2 -lpthread

CC=g++
STDLIB=-std=c++14
CFLAGS=-Wall $(STDLIB) -D_GLIBCXX_USE_CXX11_ABI=$(ABI) -O3
SOURCES=BED.cpp Expression.cpp GTF.cpp RNASeQC.cpp Metrics.cpp Fasta.cpp BamReader.cpp
SRCDIR=src
OBJECTS=$(SOURCES:.cpp=.o)
SEQFLAGS=$(STDLIB) -D_GLIBCXX_USE_CXX11_ABI=$(ABI)
SHELL=/bin/bash

rnaseqc: $(foreach file,$(OBJECTS),$(SRCDIR)/$(file)) SeqLib/lib/libseqlib.a SeqLib/lib/libhts.a
	$(CC) -O3 $(LIBRARY_PATHS) -o $@ $^ $(STATIC_LIBS) $(LIBS)

.PHONY: lib

lib: $(foreach file,$(OBJECTS),$(SRCDIR)/$(file))
	ar -rcs rnaseqc.a $^

%.o: %.cpp
	$(CC) $(CFLAGS) -I. $(INCLUDE_DIRS) -c -o $@ $<

SeqLib/lib/libseqlib.a SeqLib/lib/libhts.a:
	cd SeqLib && ./configure && make CXXFLAGS="$(SEQFLAGS)" && make install

.PHONY: clean

clean:
	rm $(wildcard $(SRCDIR)/*.o)

# The rest of the makefile consists of test cases. Run "make test" to perform all tests

.PHONY: test

test: test-version test-single test-chr1 test-downsampled test-legacy test-crams test-expected-failures
	echo Tests Complete

.PHONY: test-version

test-version: rnaseqc
	[ ! -z "$(shell ./rnaseqc --version)" ]

.PHONY: test-single

test-single: rnaseqc
	./rnaseqc test_data/single_pair.gtf test_data/single_pair.bam .test_output
	python3 test_data/approx_diff.py .test_output/single_pair.bam.metrics.tsv test_data/single_pair.output/single_pair.bam.metrics.tsv -m metrics -c single_pair.bam single_pair.bam_
	python3 test_data/approx_diff.py .test_output/single_pair.bam.gene_reads.gct <(gzcat test_data/single_pair.output/single_pair.bam.gene_reads.gct.gz) -m tables -c Counts Counts_
	python3 test_data/approx_diff.py .test_output/single_pair.bam.gene_tpm.gct <(gzcat test_data/single_pair.output/single_pair.bam.gene_tpm.gct.gz) -m tables -c TPM TPM_
	python3 test_data/approx_diff.py .test_output/single_pair.bam.exon_reads.gct <(gzcat test_data/single_pair.output/single_pair.bam.exon_reads.gct.gz) -m tables -c Counts Counts_
	python3 test_data/approx_diff.py .test_output/single_pair.bam.gene_fragments.gct <(gzcat test_data/single_pair.output/single_pair.bam.gene_fragments.gct.gz) -m tables -c Fragments Fragments_
	rm -rf .test_output

.PHONY: test-chr1

test-chr1: rnaseqc
	./rnaseqc test_data/chr1.gtf test_data/chr1.bam .test_output --coverage
	python3 test_data/approx_diff.py .test_output/chr1.bam.metrics.tsv test_data/chr1.output/chr1.bam.metrics.tsv -m metrics -c chr1.bam chr1.bam_
	python3 test_data/approx_diff.py .test_output/chr1.bam.gene_reads.gct <(gzcat test_data/chr1.output/chr1.bam.gene_reads.gct.gz) -m tables -c Counts Counts_
	python3 test_data/approx_diff.py .test_output/chr1.bam.gene_tpm.gct <(gzcat test_data/chr1.output/chr1.bam.gene_tpm.gct.gz) -m tables -c TPM TPM_
	python3 test_data/approx_diff.py .test_output/chr1.bam.exon_reads.gct <(gzcat test_data/chr1.output/chr1.bam.exon_reads.gct.gz) -m tables -c Counts Counts_
	python3 test_data/approx_diff.py .test_output/chr1.bam.gene_fragments.gct <(gzcat test_data/chr1.output/chr1.bam.gene_fragments.gct.gz) -m tables -c Fragments Fragments_
	sed s/-nan/nan/g .test_output/chr1.bam.coverage.tsv > .test_output/coverage.tsv
	python3 test_data/approx_diff.py .test_output/coverage.tsv test_data/chr1.output/chr1.bam.coverage.tsv -m metrics -c coverage_mean coverage_mean_
	python3 test_data/approx_diff.py .test_output/coverage.tsv test_data/chr1.output/chr1.bam.coverage.tsv -m metrics -c coverage_std coverage_std_
	python3 test_data/approx_diff.py .test_output/coverage.tsv test_data/chr1.output/chr1.bam.coverage.tsv -m metrics -c coverage_CV coverage_CV_
	rm -rf .test_output

.PHONY: test-downsampled

test-downsampled: rnaseqc
	./rnaseqc test_data/downsampled.gtf test_data/downsampled.bam --bed test_data/downsampled.bed --coverage .test_output
	python3 test_data/approx_diff.py .test_output/downsampled.bam.metrics.tsv test_data/downsampled.output/downsampled.bam.metrics.tsv -m metrics -c downsampled.bam downsampled.bam_
	python3 test_data/approx_diff.py .test_output/downsampled.bam.gene_reads.gct <(gzcat test_data/downsampled.output/downsampled.bam.gene_reads.gct.gz) -m tables -c Counts Counts_
	python3 test_data/approx_diff.py .test_output/downsampled.bam.gene_tpm.gct <(gzcat test_data/downsampled.output/downsampled.bam.gene_tpm.gct.gz) -m tables -c TPM TPM_
	python3 test_data/approx_diff.py .test_output/downsampled.bam.exon_reads.gct <(gzcat test_data/downsampled.output/downsampled.bam.exon_reads.gct.gz) -m tables -c Counts Counts_
	python3 test_data/approx_diff.py .test_output/downsampled.bam.gene_fragments.gct <(gzcat test_data/downsampled.output/downsampled.bam.gene_fragments.gct.gz) -m tables -c Fragments Fragments_
	sed s/-nan/nan/g .test_output/downsampled.bam.coverage.tsv > .test_output/coverage.tsv
	python3 test_data/approx_diff.py .test_output/coverage.tsv test_data/downsampled.output/downsampled.bam.coverage.tsv -m metrics -c coverage_mean coverage_mean_
	python3 test_data/approx_diff.py .test_output/coverage.tsv test_data/downsampled.output/downsampled.bam.coverage.tsv -m metrics -c coverage_std coverage_std_
	python3 test_data/approx_diff.py .test_output/coverage.tsv test_data/downsampled.output/downsampled.bam.coverage.tsv -m metrics -c coverage_CV coverage_CV_
	python3 test_data/approx_diff.py .test_output/downsampled.bam.fragmentSizes.txt test_data/downsampled.output/downsampled.bam.fragmentSizes.txt -m fragments -c Count Count_
	rm -rf .test_output

.PHONY: test-legacy

test-legacy: rnaseqc
	./rnaseqc test_data/downsampled.gtf test_data/downsampled.bam --bed test_data/downsampled.bed --coverage .test_output --legacy
	python3 test_data/approx_diff.py .test_output/downsampled.bam.metrics.tsv test_data/legacy.output/downsampled.bam.metrics.tsv -m metrics -c downsampled.bam downsampled.bam_
	python3 test_data/approx_diff.py .test_output/downsampled.bam.gene_reads.gct <(gzcat test_data/legacy.output/downsampled.bam.gene_reads.gct.gz) -m tables -c Counts Counts_
	python3 test_data/approx_diff.py .test_output/downsampled.bam.gene_tpm.gct <(gzcat test_data/legacy.output/downsampled.bam.gene_tpm.gct.gz) -m tables -c TPM TPM_
	python3 test_data/approx_diff.py .test_output/downsampled.bam.exon_reads.gct <(gzcat test_data/legacy.output/downsampled.bam.exon_reads.gct.gz) -m tables -c Counts Counts_
	python3 test_data/approx_diff.py .test_output/downsampled.bam.gene_fragments.gct <(gzcat test_data/legacy.output/downsampled.bam.gene_fragments.gct.gz) -m tables -c Fragments Fragments_
	sed s/-nan/nan/g .test_output/downsampled.bam.coverage.tsv > .test_output/coverage.tsv
	python3 test_data/approx_diff.py .test_output/coverage.tsv test_data/legacy.output/downsampled.bam.coverage.tsv -m metrics -c coverage_mean coverage_mean_
	python3 test_data/approx_diff.py .test_output/coverage.tsv test_data/legacy.output/downsampled.bam.coverage.tsv -m metrics -c coverage_std coverage_std_
	python3 test_data/approx_diff.py .test_output/coverage.tsv test_data/legacy.output/downsampled.bam.coverage.tsv -m metrics -c coverage_CV coverage_CV_
	python3 test_data/approx_diff.py .test_output/downsampled.bam.fragmentSizes.txt test_data/legacy.output/downsampled.bam.fragmentSizes.txt -m fragments -c Count Count_
	python3 test_data/approx_diff.py .test_output/downsampled.bam.gene_reads.gct <(gzcat test_data/legacy.output/legacy.gene_reads.gct.gz) -m tables -c Counts RNA-SeQC
	python3 python/rnaseqc/legacy_exon_remap.py .test_output/downsampled.bam.exon_reads.gct test_data/downsampled.gtf > /dev/null
	python3 test_data/approx_diff.py .test_output/downsampled.bam.exon_reads.gct <(gzcat test_data/legacy.output/legacy.exon_reads.gct.gz) -m tables -c Counts RNA-SeQC -t
	rm -rf .test_output

.PHONY: test-crams

test-crams: rnaseqc
	touch test_data/chr1.fasta.fai
	./rnaseqc test_data/chr1.gtf test_data/chr1.cram .test_output --coverage --fasta test_data/chr1.fasta
	python3 test_data/approx_diff.py .test_output/chr1.cram.metrics.tsv test_data/chr1.output/chr1.cram.metrics.tsv -m metrics -c chr1.cram chr1.cram_
	python3 test_data/approx_diff.py .test_output/chr1.cram.gene_reads.gct <(gzcat test_data/chr1.output/chr1.bam.gene_reads.gct.gz) -m tables -c Counts Counts_
	python3 test_data/approx_diff.py .test_output/chr1.cram.gene_tpm.gct <(gzcat test_data/chr1.output/chr1.bam.gene_tpm.gct.gz) -m tables -c TPM TPM_
	python3 test_data/approx_diff.py .test_output/chr1.cram.exon_reads.gct <(gzcat test_data/chr1.output/chr1.bam.exon_reads.gct.gz) -m tables -c Counts Counts_
	python3 test_data/approx_diff.py .test_output/chr1.cram.gene_fragments.gct <(gzcat test_data/chr1.output/chr1.bam.gene_fragments.gct.gz) -m tables -c Fragments Fragments_
	python3 test_data/approx_diff.py .test_output/chr1.cram.gc_content.tsv test_data/chr1.output/chr1.cram.gc_content.tsv -m metrics -c Count Count_
	sed s/-nan/nan/g .test_output/chr1.cram.coverage.tsv > .test_output/coverage.tsv
	python3 test_data/approx_diff.py .test_output/coverage.tsv test_data/chr1.output/chr1.bam.coverage.tsv -m metrics -c coverage_mean coverage_mean_
	python3 test_data/approx_diff.py .test_output/coverage.tsv test_data/chr1.output/chr1.bam.coverage.tsv -m metrics -c coverage_std coverage_std_
	python3 test_data/approx_diff.py .test_output/coverage.tsv test_data/chr1.output/chr1.bam.coverage.tsv -m metrics -c coverage_CV coverage_CV_
	rm -rf .test_output

.PHONY: test-expected-failures

test-expected-failures: rnaseqc
	./rnaseqc test_data/gencode.v26.collapsed.gtf test_data/downsampled.bam .test_output 2>/dev/null; test $$? -eq 11
	rm -rf .test_output
