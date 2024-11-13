#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_main as config
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import mapper
from ma_mapper import custom_cmap
from ma_mapper import sequence_alignment
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
maf_dir = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/241genomes/241-mammalian-2020v2b.maf'
e_value_table = f'{config.e_value_folder}/{subfamily}.txt'

#%%
# IMPORTANT:
# count_arg is now set to count_arg='common' by default this one counts the common base and then finds the percentage of it against total alignment
# human_ref = counts same the base as human reference then find a percentage
# coverage/total_raw = counts numbers of alignment
# common_raw = finds a raw count of the common base
# if older results need to be reproduced, use:
#maf=mapper.map_and_overlay(alignment_file, coord_file, maf_dir, data_format='maf', separated_maf = True, age_arg='calibrate', species_age_table_file=species_age_table, count_arg='human_ref')  
#%%
te_age_internal_id['meta_id'] = subfamily + '_' + te_age_internal_id.index.astype(str)
#%%
from ma_mapper import extract_maf
import importlib
importlib.reload(extract_maf)
importlib.reload(mapper)
maf=mapper.map_and_overlay(alignment_file, coord_file, maf_dir, data_format='maf', separated_maf = True, e_value_table = e_value_table, internal_id_table = te_age_internal_id, count_arg='common', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
mean_maf=mapper.normalise(alignment=alignment_filtered, mapped_data=maf)
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [maf], 
    alignment=alignment_filtered,
    show_alignment=False, 
    heatmap_color=['viridis'],
    heatmap_mode='overlay', 
    vlim = [[0,1.0]],
    show_alignment_colbar=False,
    opacity = 1.0, 
    hm_interpolation = None,
    aggregated=False,
    colorbar=True,
    colorbar_steps=[0.1], 
    aggregated_data=[mean_maf], 
    agg_colset=['green'],
    agg_ylim=[[0,1]],
    #agg_titles=['NFkB-p65 motif'], 
    #agg_ylabel=['perc_coverage'],
    hm_plot_title =f'Common base frequency\nfrom Zoonomia on THE1C MSA',
    hm_xlabel = 'position (bp)',
    hm_ylabel = 'sequences',
)
#%%
from ma_mapper import extract_maf
import importlib
importlib.reload(extract_maf)
importlib.reload(mapper)
maf=mapper.map_and_overlay(alignment_file, coord_file, maf_dir, data_format='maf', separated_maf = True, species_list=['Homo_sapiens', 'Gorilla_gorilla'], count_arg='base_count', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
#%%
import compress_pickle
output_filepath='/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/THE1C_447/the1c_maf.p'
compress_pickle.dump(maf, output_filepath, compression="lzma")
#%%
maf=compress_pickle.load(output_filepath, compression="lzma")

# %%
from ma_mapper import plots
plots.plot_experimental(data = [maf], alignment= alignment_filtered,h_cmap=['viridis'], vlim = [[0.99,1]], opacity = 1, colorbar=True, cbar_steps = [0.001])

# %%
import numpy as np
print(np.nanmin(maf),np.nanmax(maf), np.nanmean(maf), np.nanmedian(maf))
#%%
from ma_mapper import extract_maf
import importlib
importlib.reload(extract_maf)
test=extract_maf.maf_io(metadata=coord_file, maf=maf_dir,separated_maf=True, count_arg='common', e_value_table = e_value_table, internal_id_table = internal_id_sort, custom_id = True)

# %%
