#%%
import pandas as pd
import numpy as np
import sys
import os
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/stable')
import config_hg38_repeatlib as config
#%%
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
#%%
#subfamily = 'L1PA7'
#config.subfam_tally.index
for subfamily in config.subfamily_list:
    subfamily_filename = subfamily.replace('/','_')
    output_folder = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/age_redo_bed'
    output_filepath = f'{output_folder}/{subfamily_filename}.txt'
    print(f'process: {subfamily}')
    #subfamily = 'MER5C'
    #coord_internal_id_folder = config.coord_internal_id_folder
    #te_coord_tbl = f'{coord_internal_id_folder}/{subfamily_filename}.txt'
    #te_coord_df = pd.read_csv(te_coord_tbl, sep='\t')
    filtered_table = config.filtered_table
    internal_id_folder = config.internal_id_folder
    internal_id_tbl = f'{internal_id_folder}/{subfamily_filename}.txt'
    internal_id_df = pd.read_csv(internal_id_tbl, sep='\t', index_col= 0)
    te_age_folder = config.te_age_folder
    te_age_tbl = f'{te_age_folder}/{subfamily_filename}.txt'
    te_age_df = pd.read_csv(te_age_tbl, sep='\t')

    subfam_df = filtered_table[filtered_table.repName == subfamily]
    subfam_w_internal_id=pd.merge(subfam_df, internal_id_df, left_index = True, right_on = 'rmsk_index')
    subfam_w_internal_id_w_age =pd.merge(subfam_w_internal_id, te_age_df, on = 'internal_id')
    subfam_age_evi = subfam_w_internal_id_w_age[~subfam_w_internal_id_w_age.te_age.isna()]
    age_np_evi=subfam_age_evi.te_age.to_numpy()

    age_canon = np.unique(age_np_evi)
    q1, q3 = np.percentile(age_np_evi, [25, 75])
    iqr = q3 - q1
    upper_bound = q3 + (1.5 * iqr)
    adjusted_bound=find_nearest(age_canon, upper_bound)
    e_value_folder = config.e_value_folder
    e_tbl = f'{e_value_folder}/{subfamily_filename}.txt'
    e_tbl_df = pd.read_csv(e_tbl, sep='\t', low_memory=False)
    unique_id = e_tbl_df.internal_id.unique()

    e_tbl_filtered = e_tbl_df[
        ~((e_tbl_df.chr_code != '0.0') & (e_tbl_df.divergence > adjusted_bound) &
        (e_tbl_df.groupby(['internal_id', 'divergence'])['divergence'].transform('size') <= 1))
    ]
    output_folder = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/e_evi_filtered'
    output_filepath = f'{output_folder}/{subfamily}.txt'
    e_tbl_filtered.to_csv(output_filepath, sep='\t', index = False)
    print(f'{subfamily} done')
#%%
    