#%%
import sys
from turtle import width
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_main as config
from ma_mapper import mapper
from ma_mapper import custom_cmap
import numpy as np
#%%
subfamily='MER11'
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
#%%
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file,col_threshold = 0.0, col_content_threshold = 0.0, row_threshold = 0.0)
age_table_list = [f'{config.te_age_folder}/MER11A.txt',
             f'{config.te_age_folder}/MER11B.txt',
             f'{config.te_age_folder}/MER11C.txt']
age_df_list = []
for age_tbl in age_table_list:
    age_df_list.append(pd.read_csv(age_tbl, sep='\t'))
age_df=pd.concat(age_df_list)
internal_id_tbl_list = [f'{config.internal_id_folder}/MER11A.txt',
             f'{config.internal_id_folder}/MER11B.txt',
             f'{config.internal_id_folder}/MER11C.txt']
internal_id_df_list = []
for internal_id_tbl in internal_id_tbl_list:
    internal_id_df_list.append(pd.read_csv(internal_id_tbl, sep='\t').sort_values('rmsk_index'))
internal_id_sort=pd.concat(internal_id_df_list)
te_age_internal_id=internal_id_sort.merge(age_df, on='internal_id', how='left')
#%%
age_default_id = pd.DataFrame()
age_default_id['internal_id'] = subfamily + '_' + te_age_internal_id.index.astype(str)
age_default_id['group'] = te_age_internal_id.internal_id
age_default_id['te_age'] = te_age_internal_id['te_age']
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered, age_table= age_default_id)
#%%
coord_file=mapper.extract_metadata_from_alignment(alignment_file)
# %%
# %% 
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    show_alignment=True,
    alignment=alignment_filtered, 
    alignment_col='dna', 
    show_alignment_colbar=True,
    colorbar=True,
    hm_plot_title = 'MER11 MSA (filtered)',
    hm_xlabel = 'position (bp)',
    hm_ylabel = 'sequences',

    )
#%%

# %%
