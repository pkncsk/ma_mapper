#%%
import pandas as pd
import numpy as np
import os
import sys
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/stable')
import config_hg38_repeatlib as config
from concurrent.futures import ProcessPoolExecutor
#%% mm39 mouse genome setting
#%%
#subfamily = 'MER5C'
def extract_millDiv(subfamily):
    output_folder = config.kimura_distance_folder
    filtered_table = config.filtered_table
    subfamily_filename = subfamily.replace('/','_')
    print(f'process: {subfamily}')
    output_filepath = f'{output_folder}/{subfamily_filename}.txt'
    internal_id_folder = config.internal_id_folder
    internal_id_tbl = f'{internal_id_folder}/{subfamily_filename}.txt'
    internal_id_df = pd.read_csv(internal_id_tbl, sep='\t', index_col= 0)
    intersect_df = internal_id_df.merge(filtered_table, left_on='rmsk_index', right_index=True)
    def calculate_div(group):
            frag_length = group['genoEnd'] - group['genoStart']
            frag_div = frag_length * group['milliDiv']
            return frag_div.sum() / frag_length.sum()
    div_tbl = intersect_df.groupby('internal_id').apply(calculate_div,include_groups=False).reset_index(name='te_div')
    div_tbl.to_csv(output_filepath, sep='\t', index = False)
    print(f'done: {subfamily}')
#%%
# %%
def main():
    subfamily_list = config.subfamily_list
    with ProcessPoolExecutor(max_workers=40) as executor:
        results = executor.map(extract_millDiv,subfamily_list)
#%%
if __name__ == '__main__':
    main()
#%%