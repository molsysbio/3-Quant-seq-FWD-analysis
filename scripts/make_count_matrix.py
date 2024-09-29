import glob
import os
import pandas as pd
import re


# --- IMPORT FILES AND MERGE --------------------------------------------------
# get a list of the samples filepaths
# filepaths = glob.glob('subread/*_gene_level.txt')
filepaths = snakemake.input

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
# from data structure we get '/data/bla/alignment/...' in beginning --> remove
# matrix.columns = matrix.columns.str.replace('alignment/bam_dedup/', '')
matrix.rename(columns={
    col: col.split('/')[-1] for col in matrix.columns if col != 'Geneid'}, 
    inplace=True
)
# from naming of bcl2fastq we get '_S**_R1_001_dedup.bam" at the end --> remove
delimiter_pattern = r'_S\d{1,}'
matrix.rename(columns={
    col: re.split(delimiter_pattern, col)[0] for col in matrix.columns if col != 'Geneid'}, 
    inplace=True
)

# drop duplicate row entries based geneid (keep the one with highest counts)
# and de-select annotations columns
matrix['row_sum'] = matrix.sum(axis=1, numeric_only=True)
matrix = (matrix
 .sort_values('row_sum', ascending=False)
 .drop_duplicates(subset=['Geneid'])
 .sort_index()
 .drop(columns = ['row_sum'])
)

# put gene ID first two column
first_cols = ['Geneid', ]
# Get the remaining columns
remaining_columns = [col for col in matrix.columns if col not in first_cols]
# Reorder the columns
new_order = first_cols + remaining_columns
matrix = matrix[new_order]


# ---- save the count matrix --------------------------------------------------
# save count matrix as tab-delimited file
matrix.to_csv(snakemake.output.count_matrix, index=False, sep='\t')


# ---- calculate library sizes and then save them -----------------------------
(matrix
 .select_dtypes(include='number')
 .sum(axis=0)       # use axis=0 here, as you want so sum along the rows
 .sort_values()     # sorts smallest values first
 .to_csv(snakemake.output.library_sizes, sep='\t', index=True,
    header=['library_size'], index_label='sample')
)
# use index=True, as .sum returns a series, where column names are the index
# !! sorts smallest values first, so that you can see which samples are unusable !!
