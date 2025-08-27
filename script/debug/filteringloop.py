from ma_mapper import mapper
from ma_mapper import plots
from ma_mapper import custom_cmap
import numpy as np
import pandas as pd
#%%
subfamily='THE1C'
alignment_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/alignment/THE1C.fasta.aligned'
genomewide_data_filepath  = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/UCSC_phyloP_track/hg38.phyloP447way.bw'
#%%
alignment_matrix,alignment_coordinate,filter  = mapper.parse_and_filter(alignment_file=alignment_filepath,col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10, preprocess_out=True)
#%%
col_threshold = 0.10
row_threshold = 0.10
num_rows, num_cols = alignment_matrix.shape

    # Start with all rows and columns included
current_row_indices = np.arange(num_rows)
current_col_indices = np.arange(num_cols)

while True:
    submatrix = alignment_matrix[np.ix_(current_row_indices, current_col_indices)]

    # Compute column nonzero ratios
    col_nonzero_counts = np.count_nonzero(submatrix, axis=0)
    col_nonzero_ratio = col_nonzero_counts / submatrix.shape[0]

    # Filter columns by threshold
    cols_to_keep_mask = col_nonzero_ratio >= col_threshold
    new_col_indices = current_col_indices[cols_to_keep_mask]

    # Compute row nonzero ratios (on filtered columns)
    submatrix_rows_filtered = alignment_matrix[np.ix_(current_row_indices, new_col_indices)]
    row_nonzero_counts = np.count_nonzero(submatrix_rows_filtered, axis=1)
    row_nonzero_ratio = row_nonzero_counts / submatrix_rows_filtered.shape[1]

    # Filter rows by threshold
    rows_to_keep_mask = row_nonzero_ratio >= row_threshold
    new_row_indices = current_row_indices[rows_to_keep_mask]

    # Check for convergence: if indices didn’t change, break
    if (np.array_equal(new_row_indices, current_row_indices) and
        np.array_equal(new_col_indices, current_col_indices)):
        break

    # Update indices for next iteration
    current_row_indices = new_row_indices
    current_col_indices = new_col_indices
# %%
