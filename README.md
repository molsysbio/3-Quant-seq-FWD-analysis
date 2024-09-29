
# QuantSeq FWD-UMI Data Analysis Pipeline

This snakemake-based pipeline is specific for QuantSeq FWD libraries that contain Unique Molecular Identifiers (UMIs) in Read 1. It is designed to be used on the BIH cluster. Starting from fastq files, reads are aligned using STAR, deduplicated based on UMI and counted by subread's featureCounts. Quality reports of the unprocessed and processed fastq files are generated with fastqc. Finally, a count matrix is generated. 

To get started, you only need to provide a file path where all the fastq files from the 3' Quant-seq experiment are stored. The pipeline then automatically selects all forward read files by filtering for `*_R1_001.fastq.gz`. This catches all relevant files if you demultiplexed your data with bcl2fastq. If you renamed the fastq files, you can change the filtering behaviour in line 20 in the Snakefile.

The folder `example_fastq_quant-seq` contains three fastq files with only 5000 reads from a 3' Quant-seq FWD UMI experiment. To start a test run, the filepaths can be left as is. Running the whole pipeline on the example data should take only ~5 minutes. Installing the conda enviroments can take much longer though.

## Prerequisites

### 1. Conda and Snakemake Installation
The pipeline relies on Conda for managing dependencies and Snakemake for workflow execution. To install Snakemake, please make sure you have Conda installed, then run:

```bash
conda install -c bioconda snakemake
```

**Important**: Ensure that the Conda channel priority is set to **flexible** (default). You can set this by running:

```bash
conda config --set channel_priority flexible
```

### 2. STAR Precomputed Genomes

The pipeline requires precomputed STAR genome indices. These should be provided in the configuration file (`config.yaml`). To generate a STAR genome index, you need a genome fasta file and the corresponding annotation GTF file (e.g. provided by [GENCODE](https://www.gencodegenes.org/human/release_46.html)). The following command generates a STAR index. Make sure the STAR version is *2.7.1a*.

```bash
STAR --runMode genomeGenerate --genomeDir /path/to/STAR_INDEX/ --genomeFastaFiles /path/to/genome.fa --sjdbGTFfile /path/to/annotation.gtf --sjdbOverhang 75 --runThreadN 16 --limitGenomeGenerateRAM 100000000000
```

### 3. UMI Barcode Pattern

The pipeline is designed for data where UMIs are found in Read 1, and the UMI extraction rule is based on a barcode pattern. This can be specified in the `config.yaml` file. For example, a 6 bp UMI + 4 bp spacer would be represented as 10 'N's:

```yaml
bc_pattern: NNNNNNNNNN
```

Modify this based on your actual UMI and spacer configuration. The default of 6 bp UMI + 4 bp spacer is very common though.

## Configuration

Before running the pipeline, you must adjust the `config.yaml` file to specify the paths to your data and other parameters. Here are the key settings in the file:

- **fastq_dir**: Path to the directory containing FASTQ files (with a trailing slash).
- **count_matrix**: Name of the output file for the count matrix.
- **library_sizes**: Name of the output file for the summary of library sizes.
- **tmp**: Path for temporary files, set to a scratch directory.
- **subread_features**: Feature level for gene quantification (e.g., `gene`, `transcript`, `exon`).

```yaml
fastq_dir: /path/to/your/fastq/files/
count_matrix: count_matrix.tsv
library_sizes: library_sizes.tsv
tmp: "$HOME/scratch/"
subread_features: gene
```

Ensure that you customize these fields before running the pipeline.

## Running the Pipeline

Once all dependencies are installed and the configuration file is set, you can run the pipeline with the following command:

```bash
snakemake --profile helper_files/brecht_profile/ -j 30
```

Replace 30 with the number of jobs you want to run in parallel. 

For dry-run mode (to verify that everything is set up correctly), you can use:

```bash
snakemake -n
```

## Output Files

- **Count Matrix**: A tab-delimited file containing gene counts across samples (by default `count_matrix.tsv`).
- **Library Sizes**: A summary of library sizes across samples (by default `library_sizes.tsv`).

## Notes

- Ensure that all paths provided in the `config.yaml` file are correct and accessible from your computing environment.
- Adjust UMI extraction patterns and STAR genome settings according to your specific dataset.
