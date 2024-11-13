#%%
import pandas as pd
import os
import sys
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
from concurrent.futures import ProcessPoolExecutor
#%%
#%%def make_te_age_bed(subfamily):
subfamily = 'THE1C'
subfamily_filename = subfamily.replace('/','_')
print(f'process: {subfamily}')
#subfamily = 'MER5C'
#coord_internal_id_folder = config.coord_internal_id_folder
#te_coord_tbl = f'{coord_internal_id_folder}/{subfamily_filename}.txt'
#te_coord_df = pd.read_csv(te_coord_tbl, sep='\t')
internal_id_folder = config.internal_id_folder
internal_id_tbl = f'{internal_id_folder}/{subfamily_filename}.txt'
internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
te_age_folder = config.te_age_folder
te_age_tbl = f'{te_age_folder}/{subfamily_filename}.txt'
te_age_df = pd.read_csv(te_age_tbl, sep='\t')
te_div_folder = config.kimura_distance_folder
te_div_tbl = f'{te_div_folder}/{subfamily_filename}.txt'
te_div_df = pd.read_csv(te_div_tbl, sep='\t')
#%%
internal_id_div = internal_id_df.merge(te_div_df, on = 'internal_id')
internal_id_age=internal_id_div.merge(te_age_df, on = 'internal_id')
# %%
