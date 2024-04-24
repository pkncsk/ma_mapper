#%%
import sys
import pandas as pd
import numpy as np
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import mapper
from ma_mapper._unused import fetch_data
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import gcf
import logging
logging.getLogger('matplotlib').setLevel(logging.WARNING)
#%%
input_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11.fasta.aligned'
aligned_parsed = mapper.parse_alignment(input_filepath, save_to_file= False)
metadata_aligned = mapper.extract_metadata_from_alignment(input_filepath)
# %%
metadata_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11_coord_with_id_age.txt'
#%%
metadata_df = pd.read_csv(metadata_filepath, sep='\t')
original_order = metadata_df.iloc[:,4].unique()
#%% make new metadata 
age_table_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/species241_info.tsv'

maf_mapped=fetch_data.fetch_maf(metadata_input= metadata_filepath, maf_input='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/241genomes/241-mammalian-2020v2b.maf', separated_maf = True,target_species = 'Homo_sapiens', custom_id= True, count_arg= 'common',age_arg = 'calibrate', age_table_file= age_table_filepath)
#%%
maf_mapped_sorted = []
for idx, row in metadata_aligned.iterrows():
    maf_mapped_sorted.append(maf_mapped[np.where(original_order == row.id)[0][0]])
metadata_df_sorted=metadata_aligned.merge(metadata_df[['id','te_age']].drop_duplicates
(), on='id', how='left')
#%%
filters=mapper.create_filter(aligned_parsed)
row_filter = filters[0]
col_filter = filters[1]
#%%
aligned_maf_overlay=mapper.map_data(maf_mapped_sorted, aligned_parsed, filters = filters)
metadata_df_sorted_filtered=metadata_df_sorted.iloc[row_filter,:]
#%%
plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600
fig = plt.figure(figsize=(6,5))
grid = fig.add_gridspec(nrows = 1, ncols = 1, hspace=0)
ax0 = fig.add_subplot(grid[0,0])
heatmap_core = ax0.imshow(aligned_maf_overlay, aspect = 'auto', cmap='viridis', vmax = 1)

cbar = ax0.figure.colorbar(heatmap_core, ax=ax0)
plt.show()
#%%
aligned_maf_overlay
# %%
