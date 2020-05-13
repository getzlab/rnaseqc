# RNA-SeQC Python utilities

This module contains utility code for RNA-SeQC

## Installing

* From pip: `pip install rnaseqc`
* From the git repo: `pip install -e python` (Invoke from root of git repo)

## Usage

This does not install a console entrypoint. You can invoke the utilities in one of three ways:

* From the main module: `python3 -m rnaseqc ...`
* Calling the target module: `python3 -m rnaseqc.example ...`
* Calling scripts directly: `python3 python/rnaseqc/example.py`

## Utilities

The `rnaseqc` module contains 5 main utilities. To get more help with each utility,
invoke the utility with the `-h` or `--help` option

### Aggregation

Aggregates RNA-SeQC outputs from multiple samples

```
python3 -m rnaseqc aggregate [-h] [--parquet] [-o OUTPUT_DIR] results_dir prefix
```

### Jupyter Notebooks

Creates a jupyter notebook with several figures for comparing samples

```
python3 -m rnaseqc notebook [-h] [-t TPM] [-i INSERT_SIZE] [-c COHORT] [-d DATE] metrics output
```

### Figures

Generates figures from an aggregated RNA-SeQC metrics table

```
python3 -m rnaseqc report [-h] [--tpm TPM] [--insert-size INSERT_SIZE] [--cohort COHORT] [--output-dir OUTPUT_DIR] [--dpi DPI] metrics prefix
```

### Insert Size distributions

Generates a BED file with intervals used by RNA-SeQC for estimating a sample's insert size distribution

```
python3 -m rnaseqc insert-size [-h] [--min-length MIN_LENGTH] [--min-mappability MIN_MAPPABILITY] [--output-dir OUTPUT_DIR] gtf_path mappability_bigwig prefix
```

### Exon remapping

Convert exon names in an `*.exon_reads.gct` file from RNA-SeQC 2.X.X to match names
as reported by RNA-SeQC 1.1.9

```
python3 -m rnaseqc legacy-exons gct gtf
```
