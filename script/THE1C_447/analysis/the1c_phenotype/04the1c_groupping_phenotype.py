#%%
import sys
from turtle import width
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
from ma_mapper import mapper
from ma_mapper import custom_cmap
import numpy as np
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
#%%
coord_file=mapper.extract_metadata_from_alignment(alignment_file)
#%%
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/zoonomia447/hg38.phyloP447way.bw'
phylop_447=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', custom_id=True)
# %%
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/AP-1(bZIP).bed'
ap1=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True)
#%%
mean_phylop=mapper.normalise(alignment=alignment_filtered, mapped_data=phylop_447)
mean_ap1=mapper.normalise(alignment=alignment_filtered, mapped_data=ap1)
# %%
#find peaks
import scipy
peaks, _ = scipy.signal.find_peaks(mean_ap1, width = 8)
#highest_peak_index = peaks[np.argmax(mean_znf808[peaks])]
binding_indices = np.unique(np.where(ap1[:, peaks] != 0)[0])
nonbinding_indices=list(set(np.arange(ap1.shape[0])) - set(binding_indices))
#%%
# %%
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
phylop_sorted = np.vstack((phylop_bind,phylop_nonbind))
ap1_bind = ap1[binding_indices]
ap1_nonbind = ap1[nonbinding_indices]
ap1_sorted = np.vstack((ap1_bind,ap1_nonbind))
te_age_sorted=metadata_age.iloc[np.concatenate((binding_indices,nonbinding_indices))].te_age.fillna(0)
# %%
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
stat_v, p_value = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#%%
anno_label=metadata_age.fillna(0).te_age.sort_values().unique()
# %% 
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(data = [phylop_sorted,ap1_sorted,], alignment=alignment_filtered, h_cmap=[custom_cmap.vlag_r_mpl,'Greens'], vlim = [[-0.5,0.5],[0,0.0005],[0,0.0005]], opacity = 0.5, annotation=True, anno_col = ['Blues'],heatmap_mode='spread_horizontal', annotation_data=[te_age_sorted], aggregated=True, aggregated_data=[mean_phylop,-np.log(p_value)], a_colset=['grey','grey'],agg_ylim=[[-1,1],[0,50]], colorbar=True,cbar_steps = [0.1,0.0001], anno_cbar_label=[anno_label], agg_titles=['mean_phyloP','log p-value H1: bind greater\nthan nonbind'], anno_title=['age'],anno_cbar_title=['MYA'])
#%%