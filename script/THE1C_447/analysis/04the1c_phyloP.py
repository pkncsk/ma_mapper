#%%
from configparser import Interpolation
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
subfamily = 'THE1C'
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file, col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
age_table = f'{config.te_age_folder}/{subfamily}.txt'
#%%
age_df = pd.read_csv(age_table, sep='\t')
internal_id_tbl = f'{config.internal_id_folder}/{subfamily}.txt'
internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
internal_id_sort = internal_id_df.sort_values('rmsk_index')
te_age_internal_id=internal_id_sort.merge(age_df, on='internal_id', how='left')
#%%
age_default_id = pd.DataFrame()
age_default_id['internal_id'] = subfamily + '_' + te_age_internal_id.index.astype(str)
age_default_id['te_age'] = te_age_internal_id['te_age']
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered, age_table=age_default_id)
#%%
coord_file=mapper.extract_metadata_from_alignment(alignment_file)
#%%
bigwig_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_bigwig/241-mammalian-2020v2.bigWig'
phylop=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
#%%
alignment_extended, metadata_extended= mapper.parse_and_filter(alignment_file, col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10,extension_length=500, source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_fasta/hg38.fa')
#%%
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/zoonomia447/hg38.phyloP447way.bw'
phylop_447=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
# %%
#%%
mean_phylop_447=mapper.normalise(alignment=alignment_extended, mapped_data=phylop_447)
mean_phylop=mapper.normalise(alignment=alignment_filtered, mapped_data=phylop)
#%%
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6,6))
ax.scatter(mean_phylop, mean_phylop_447, alpha=1,s=2, color='black')
# Fit a line of best fit
from scipy import stats
res = stats.linregress(mean_phylop, mean_phylop_447)
# Line of best fit
ax.plot(mean_phylop, res.intercept + res.slope*mean_phylop, color='red', label='Line of Best Fit',linewidth=1, alpha=0.2)
ax.set_ylim(-1.2,0.5)
ax.set_xlim(-1.2,0.5)
ax.set_xlabel('phyloP from Zoonomia (241 species)')
ax.set_ylabel('phyloP from Cactus (447 species)')

plt.show()

# %% 
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(heatmap=True,
    data = [phylop_447,], 
    alignment=alignment_filtered, 
    heatmap_color=[custom_cmap.vlag_r_mpl],
    heatmap_mode='overlay', 
    vlim = [[-0.5,0.5],[0,0.0005],], 
    opacity = 1.0, 
    #hm_zero_col_centering = True,
    #aggregated=True, 
    #aggregated_data=[mean_phylop_447], 
    #agg_colset=['grey','grey','red'],
    #agg_ylim=[[-1.2,1.2]],
    #agg_titles=['phyloP_447',], 
    colorbar=True,
    colorbar_steps = [0.1,0.0001],  
    agg_major_tick=50,
    agg_h =10,
    hm_plot_title =f'phyloP from cactus dataset on THE1C MSA',
    hm_xlabel = 'position (bp)',
    hm_ylabel = 'sequences',
    )
#%%

# %%
