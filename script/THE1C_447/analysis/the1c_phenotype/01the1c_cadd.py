#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import extract_bigwig, mapper, sequence_alignment
from ma_mapper import custom_cmap
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
import compress_pickle
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
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/phenotype_tracks/cadd/a.bw'
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
caddA=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', custom_id=True)
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/phenotype_tracks/cadd/t.bw'
caddT=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', custom_id=True)
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/phenotype_tracks/cadd/c.bw'
caddC=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', custom_id=True)
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/phenotype_tracks/cadd/g.bw'
caddG=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', custom_id=True)
#%%
mean_caddA=mapper.normalise(alignment=alignment_filtered, mapped_data=caddA)
mean_caddT=mapper.normalise(alignment=alignment_filtered, mapped_data=caddT)
mean_caddC=mapper.normalise(alignment=alignment_filtered, mapped_data=caddC)
mean_caddG=mapper.normalise(alignment=alignment_filtered, mapped_data=caddG)
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [caddA,caddT,caddC,caddG], 
    alignment=alignment_filtered, 
    heatmap_color=['Greens','Blues','Oranges','Reds'], 
    heatmap_mode = 'spread_horizontal', 
    vlim = [[0,20]], 
    opacity = 0.5, 
    aggregated=True, 
    aggregated_data=[mean_caddA,mean_caddT,mean_caddC,mean_caddG], 
    agg_colset=['green','blue','orange','red'], 
    agg_ylim=[[0,5]])
# %%
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/AP-1(bZIP).bed'
ap1=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=False)
#%%
mean_ap1=mapper.normalise(alignment=alignment_filtered, mapped_data=ap1)
#%%
import scipy
import numpy as np
peaks, _ = scipy.signal.find_peaks(mean_ap1, width = 8)
#highest_peak_index = peaks[np.argmax(mean_znf808[peaks])]
binding_indices = np.unique(np.where(ap1[:, peaks] != 0)[0])
nonbinding_indices=list(set(np.arange(ap1.shape[0])) - set(binding_indices))
# %%
caddA_bind = caddA[binding_indices]
caddA_nonbind = caddA[nonbinding_indices]
caddA_sorted = np.vstack((caddA_bind,caddA_nonbind))
caddT_bind = caddT[binding_indices]
caddT_nonbind = caddT[nonbinding_indices]
caddT_sorted = np.vstack((caddT_bind,caddT_nonbind))
caddC_bind = caddC[binding_indices]
caddC_nonbind = caddC[nonbinding_indices]
caddC_sorted = np.vstack((caddC_bind,caddC_nonbind))
caddG_bind = caddG[binding_indices]
caddG_nonbind = caddG[nonbinding_indices]
caddG_sorted = np.vstack((caddG_bind,caddG_nonbind))
ap1_bind = ap1[binding_indices]
ap1_nonbind = ap1[nonbinding_indices]
ap1_sorted = np.vstack((ap1_bind,ap1_nonbind))
te_age_sorted=metadata_age.iloc[np.concatenate((binding_indices,nonbinding_indices))].te_age.fillna(0)
# %%
alignemnt_bind=alignment_filtered[binding_indices]
caddA_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
caddA_nonbind[alignemnt_nonbind == 0] = np.nan
caddA_stat_v_greater, caddA_p_value_greater = scipy.stats.mannwhitneyu(caddA_bind,caddA_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
caddA_stat_v_less, caddA_p_value_less = scipy.stats.mannwhitneyu(caddA_bind,caddA_nonbind, axis =0,nan_policy='omit', alternative = 'less')
# %%
alignemnt_bind=alignment_filtered[binding_indices]
caddT_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
caddT_nonbind[alignemnt_nonbind == 0] = np.nan
caddT_stat_v_greater, caddT_p_value_greater = scipy.stats.mannwhitneyu(caddT_bind,caddT_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
caddT_stat_v_less, caddT_p_value_less = scipy.stats.mannwhitneyu(caddT_bind,caddT_nonbind, axis =0,nan_policy='omit', alternative = 'less')
# %%
alignemnt_bind=alignment_filtered[binding_indices]
caddC_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
caddC_nonbind[alignemnt_nonbind == 0] = np.nan
caddC_stat_v_greater, caddC_p_value_greater = scipy.stats.mannwhitneyu(caddC_bind,caddC_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
caddC_stat_v_less, caddC_p_value_less = scipy.stats.mannwhitneyu(caddC_bind,caddC_nonbind, axis =0,nan_policy='omit', alternative = 'less')
# %%
alignemnt_bind=alignment_filtered[binding_indices]
caddG_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
caddG_nonbind[alignemnt_nonbind == 0] = np.nan
caddG_stat_v_greater, caddG_p_value_greater = scipy.stats.mannwhitneyu(caddG_bind,caddG_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
caddG_stat_v_less, caddG_p_value_less = scipy.stats.mannwhitneyu(caddG_bind,caddG_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#%%
anno_label=metadata_age.fillna(0).te_age.sort_values().unique()
# %% 
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [caddA_sorted,caddT_sorted,caddC_sorted,caddG_sorted], 
    alignment=alignment_filtered, 
    heatmap_color=['Greens','Blues','Oranges','Reds'],
    heatmap_mode='spread_horizontal', 
    vlim = [[0,20]], 
    opacity = 0.5, 
    annotation=True, 
    anno_col = ['Blues'], 
    annotation_data=[te_age_sorted],
    anno_title=['age'],
    anno_cbar_title=['MYA'], 
    anno_cbar_label=[anno_label], 
    aggregated=True, 
    agg_h=5,
    aggregated_data=[-np.log(caddA_p_value_greater),np.log(caddA_p_value_less),-np.log(caddT_p_value_greater),np.log(caddT_p_value_less),-np.log(caddC_p_value_greater),np.log(caddC_p_value_less),-np.log(caddG_p_value_greater),np.log(caddG_p_value_less)], 
    agg_colset=['blue','red','blue','red','blue','red','blue','red'], 
    agg_titles=['A logP p>a','A logP p<a','T logP p>a','T logP p<a','C logP p>a','C logP p<a','G logP p>a','G logP p<a'], 
    )
# %%
plots.plot_experimental(data = [caddA,caddT,caddC,caddG], alignment=alignment_filtered, heatmap_color=['Greens','Blues','Oranges','Reds'], heatmap_mode = 'spread_horizontal', vlim = [[0,20]], opacity = 0.5, aggregated=True, aggregated_data=[mean_caddA,mean_caddT,mean_caddC,mean_caddG], agg_colset=['green','blue','orange','red'], agg_ylim=[[0,5]])