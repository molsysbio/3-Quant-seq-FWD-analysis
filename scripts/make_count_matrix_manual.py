import glob
import os
import pandas as pd


# --- IMPORT FILES AND MERGE --------------------------------------------------
# get a list of the samples filepaths
filepaths = glob.glob('../subread/*_gene_level.txt')
# filepaths = snakemake.input

# start count matrix with first count file, to which the other samples' counts will be joined
matrix = pd.read_table(filepaths[0], sep=' ')
# sort it by gene id
matrix.sort_values(by=['Geneid'], inplace=True)

# join other files to inital sample
for file in filepaths[1:]:
    sample = pd.read_table(file, sep=' ')
    matrix = pd.merge(matrix, sample, on='Geneid', how='left', sort=False)

# change the Geneid column to exclude version number of genes
matrix.loc[:,'Geneid'] = matrix['Geneid'].str.split('.', expand=True)[0]
# !! this is necessary when aligning to Gencode genomes !!

# sort columns
matrix.sort_index(axis=1, inplace=True)

# make column names nicer
# from data structure we get 'alignment/...' in beginning --> remove
matrix.columns = matrix.columns.str.replace('alignment/bam_dedup/', '')
# from naming of bcl2fastq we get '_S**_R1_001_dedup.bam" at the end --> remove
matrix.rename(columns={
    col: '_'.join(col.split('_')[:-4]) for col in matrix.columns if col != 'Geneid'}, 
    inplace=True
)


# ---- now directly do the filtering out of rRNA -------------
# get annotations
# annotations = pd.read_csv(snakemake.params.anno)
annotations = pd.read_csv(
    '~/work/databases/annotations/annotations_hsapiens_GRCh38_p13.csv')
# join, filter out rRNA and tRNA
matrix = (matrix
 .merge(annotations, on='Geneid', how='left', sort=False)
 .query('not(\
    gene_biotype == "rRNA" | gene_biotype == "rRNA_pseudogene" | \
    gene_biotype == "processed_pseudogene" | gene_biotype == "misc_RNA" | \
    gene_biotype == "ribozyme")')
)
# drop duplicate row entries based geneid (keep the one with highest counts)
# and de-select annotations columns
matrix['row_sum'] = matrix.sum(axis=1, numeric_only=True)
matrix = (matrix
 .sort_values('row_sum', ascending=False)
 .drop_duplicates(subset=['Geneid'])
 .sort_index()
 .drop(columns = ['symbol', 'gene_biotype', 'row_sum'])
)

# save count matrix to csv
matrix.to_csv('../data/count_matrix_agne_MCF10A_mutants.tsv', index=False, sep='\t')
# matrix.to_csv(snakemake.output.count_matrix, index=False, sep='\t')


# ---- calculate library sizes and then save them -----------------
(matrix
 .select_dtypes(include='number')
 .sum(axis=0)
 .sort_values()
 .to_csv('../data/library_sizes_agne_MCF10A_mutants.tsv', sep='\t', index=True,
#  .to_csv(snakemake.output.library_sizes, sep='\t', index=True,
    header=['library_size'], index_label='sample')
)
# use axis=0 here, as you want so sum along the rows (axis=1 would give sum per row)
# use index=True, as .sum returns a series, where column names are the index
# sorts smallest values first, so that you can see which samples are unusable