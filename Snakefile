# snakemake file for 3' Quant-seq pipeline (read alignment and counting)
# all filepaths and configurations for this pipeline can be set in the config file
# ./config.yaml

# to run the Snakefile use the following command:
# "snakemake --profile cubi-v1 -j 30"
# or if you want it faster
# "snakemake --profile helper_files/brecht_profile/ -j 30"

import glob
import os


################################################################################
### preliminary part
################################################################################
configfile: "config.yaml"

# use snakemakes function to collect filenames without extension in fastq dir
(SAMPLES,) = glob_wildcards(os.path.join(config['fastq_dir'], "{id}_R1_001.fastq.gz"))

######### path variables for annotations and reference genomes
# bbduk needs reference fa files
ADAPTER_SEQUENCES = os.path.join(config['helper_files'], 'adapters_for_bbduk/')
# here, STAR genome is already precomputed (check for spdjOverhang)
STAR_X = config['star_sjdb_50']
# featureCounts needs the GTF file (same one as used for generating STAR genome)
GTF = os.path.join(config['star_sjdb_50'], 'annotation.gtf')

# expand temporary file path
TMP = os.path.expandvars(config['tmp'])


################################################################################
### Snakemake part
################################################################################
# ------ rule all --------------------------------------------------------------
# initial rule to get snakemake starting (only input given, no actual rule)
# points to the files that should be the end results
# snakemake then works its way backwards
rule all:
    input:
        expand("alignment/fast_qc_before/{sample}", sample=SAMPLES),
        expand("alignment/fast_qc_after/{sample}", sample=SAMPLES),
        config['count_matrix'],
        config['library_sizes'],


# give count matrix file as input (--> runs alignment etc),
# but also fastqc checks explicitly, since they are not used downstream


# ------ fastqc quality check number 1 -----------------------------------------
# check fastq quality before any processing steps
rule fastqc1:
    threads: 4
    resources:
        mem_mb=2000,
        time="04:00:00",
    input:
        fq=config['fastq_dir'] + "{sample}_R1_001.fastq.gz",
    output:
        directory("alignment/fast_qc_before/{sample}"),
    conda:
        "helper_files/conda_envs/fastqc.yaml"
    shell:
        "mkdir -p {output};"
        "fastqc {input} --outdir={output} --extract --threads {threads}"


# ------ adapter trimming ------------------------------------------------------
# trim polyA tails and adapters with bbduk
# ref files need to be located, but are also in the installation of bbmap
rule trimming:
    threads: 8
    resources:
        mem_mb=2000,
        time="04:00:00",
    input:
        fq=config['fastq_dir'] + "{sample}_R1_001.fastq.gz",
        fa1=ADAPTER_SEQUENCES + "polyA.fa.gz",
        fa2=ADAPTER_SEQUENCES + "truseq.fa.gz",
    output:
        os.path.join(TMP, "alignment/trimmed/{sample}_trimmed.fastq.gz"),
    conda:
        "helper_files/conda_envs/bbmap.yaml"
    shell:
        r"""
        bbduk.sh \
        in={input.fq} out={output} \
        ref={input.fa1},{input.fa2} \
        k=13 ktrim=r useshortkmers=t mink=5 \
        qtrim=r trimq=10 minlength=20 t=4
        """
# t=X option limits only number of worker threads, input and output have their own
# --> choose more threads in rule (8), but only 4 in command


# ------ extract UMIs ----------------------------------------------------------
# use umi_tools to extract UMI in fastq files
rule extract_umis:
    threads: 6
    resources:
        mem_mb=5000,
        time="04:00:00",
    input:
        os.path.join(TMP, "alignment/trimmed/{sample}_trimmed.fastq.gz"),
    output:
        os.path.join(TMP, "alignment/extracted/{sample}_extracted.fastq.gz"),
    params:
        tmp=TMP,
        bc_pattern=config['bc_pattern'],
    conda:
        "helper_files/conda_envs/umitools_v2.yaml"
    shell:
        r"""
        mkdir -p alignment/extracted/logs/

        umi_tools extract \
        --stdin={input} \
        --bc-pattern={params.bc_pattern} \
        --log=alignment/extracted/logs/{wildcards.sample}.log \
        --temp-dir={params.tmp} \
        --stdout={output}
        """
# important parameter is --bc-pattern
# here, UMIs are 6bp long + 4bp spacer (TATA) --> 10 total i.e. NNNNNNNNNN


# ------ fastqc quality check number 2 -----------------------------------------
# check fastq quality after using bbduk and umi_tools extract
rule fastqc2:
    threads: 4
    resources:
        mem_mb=2000,
        time="04:00:00",
    input:
        os.path.join(TMP, "alignment/extracted/{sample}_extracted.fastq.gz"),
    output:
        directory("alignment/fast_qc_after/{sample}"),
    conda:
        "helper_files/conda_envs/fastqc.yaml"
    shell:
        "mkdir -p {output};"
        "fastqc {input} --outdir={output} --extract --threads {threads}"


# ------ read alignment with STAR ----------------------------------------------
# align the filtered fastq files
# the fastq files from 1st and 2nd sequencing run are combined here to one BAM
rule star_align:
    threads: 16
    resources:
        mem_mb=32000,
        time="04:00:00",
    input:
        fq=os.path.join(TMP, "alignment/extracted/{sample}_extracted.fastq.gz"),
        sx=STAR_X,
    output:
        os.path.join(TMP, "alignment/bam/{sample}.bam"),
    conda:
        "helper_files/conda_envs/star.yaml"
    shell:
        r"""
        STAR \
        --genomeDir {input.sx} \
        --readFilesIn {input.fq} \
        --outFileNamePrefix alignment/bam/{wildcards.sample}_ \
        --runThreadN {threads} \
        --readFilesCommand zcat \
        --limitBAMsortRAM 30000000000 \
        --outSAMtype BAM SortedByCoordinate

        mv alignment/bam/{wildcards.sample}_Aligned.sortedByCoord.out.bam {output}
        """
# if input fastq files are compressed, use --readFilesCommand zcat
# loading options: --genomeLoad NoSharedMemory \ is standard;
# --genomeload LoadAndKeep \ exists, but many Linux servers don't allow
# a large chunck of memory to be shared --> safer to use standard option


# ------ BAM indexing ----------------------------------------------------------
# BAM indexing is required for deduplicating
rule bam_index:
    threads: 4
    resources:
        mem_mb=4000,
        time="04:00:00",
    input:
        os.path.join(TMP, "alignment/bam/{sample}.bam"),
    output:
        os.path.join(TMP, "alignment/bam/{sample}.bam.bai"),
    conda:
        "helper_files/conda_envs/samtools.yaml"
    shell:
        "samtools index {input} > {output}"


# ------ deduplication ---------------------------------------------------------
# deduplicate reads based on the umi
rule dedup_umis:
    threads: 4
    resources:
        mem_mb=8000,
        time="04:00:00",
    input:
        [os.path.join(TMP, "alignment/bam/{sample}.bam"), 
         os.path.join(TMP, "alignment/bam/{sample}.bam.bai")],
    output:
        os.path.join(TMP, "alignment/bam_dedup/{sample}_dedup.bam"),
    params:
        tmp=TMP,
    conda:
        "helper_files/conda_envs/umitools_v2.yaml"
    shell:
        r"""
        umi_tools dedup \
        --stdin={input} \
        --stdout={output} \
        --temp-dir={params.tmp}
        """


# ------ read counting ---------------------------------------------------------
# quantify the mapped reads based on features (exons) or meta-features (genes)
# subread ignores reads that overlap with multiple features (-O option not used)
rule feature_counts:
    threads: 8
    resources:
        mem_mb=4000,
        time="04:00:00",
    input:
        bam=os.path.join(TMP, "alignment/bam_dedup/{sample}_dedup.bam"),
        gtf=GTF,
    output:
        "subread/{sample}_{feature}_level.txt",
    conda:
        "helper_files/conda_envs/subread.yaml"
    shell:
        r"""
        featureCounts \
        -g gene_id \
        -T {threads} \
        -t {wildcards.feature} \
        -a {input.gtf} \
        -o subread/{wildcards.sample}_{wildcards.feature}_level_gene_id_full.txt \
        {input.bam}

        awk '{{print $1,$7}}' \
        subread/{wildcards.sample}_{wildcards.feature}_level_gene_id_full.txt | \
        tail -n +2 > {output}
        """
# -t specifies feature types,
# options: intron, exon (default), transcript, gene
# here we use gene
# output txt files: columns are separated by one whitespace i.e. " "


# ------ make the final count matrix from individual sample count data ---------
rule count_matrix:
    threads: 1
    resources:
        mem_mb=4000,
        time="04:00:00",
    input:
        expand(
            "subread/{sample}_{feature}_level.txt",
            sample=SAMPLES,
            feature=config['subread_features'],
        ),
    output:
        count_matrix=config['count_matrix'],
        library_sizes=config['library_sizes'],
    conda:
        "helper_files/conda_envs/py-env_small_versions.yaml"
    script:
        "scripts/make_count_matrix.py"
