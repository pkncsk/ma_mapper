#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
from ma_mapper import mapper
#%%
subfamily = 'THE1C'
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
age_table = f'{config.te_age_folder}/{subfamily}.txt'
age_df = pd.read_csv(age_table, sep='\t')
#%%
internal_id_tbl = f'{config.internal_id_folder}/{subfamily}.txt'
internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
internal_id_sort = internal_id_df.sort_values('rmsk_index')
te_age_internal_id=internal_id_sort.merge(age_df, on='internal_id', how='left')
#%%
age_default_id = pd.DataFrame()
age_default_id['internal_id'] = subfamily + '_' + te_age_internal_id.index.astype(str)
age_default_id['te_age'] = te_age_internal_id['te_age']
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered, age_table=age_default_id)
# %%
coord_file=mapper.extract_metadata_from_alignment(alignment_file)
#%%
bam_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/KZFP-bam_hg38/znf267.sorted.bam'
#%%
bam_forward=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_forward', custom_id=True)
bam_reverse=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_reverse', custom_id=True)
bam_min=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_min', custom_id=True)
#%%
mean_forward=mapper.normalise(alignment=alignment_filtered, mapped_data=bam_forward)
mean_reverse=mapper.normalise(alignment=alignment_filtered, mapped_data=bam_reverse)

#%%
import importlib 
from ma_mapper import plots
importlib.reload(plots)
plots.plot_experimental(data = [bam_forward,bam_reverse], h_cmap=['Blues','Reds'], vlim = [[0,0.1],[0,0.1]], hm_opacity = 1,colorbar=True, cbar_steps = [0.01,0.01], aggregated = True, a_colset = [['blue','red']], aggregated_data=[[mean_forward, mean_reverse]], ag_opacity = 0.5)
#%%
importlib.reload(plots)
plots.plot_experimental(data = [bam_forward,bam_reverse], h_cmap=['Blues','Reds'], vlim = [[0,0.1],[0,0.1]], heatmap_mode='spread_vertical' ,hm_opacity = 1,colorbar=True, cbar_steps = [0.01,0.01])
# %%
