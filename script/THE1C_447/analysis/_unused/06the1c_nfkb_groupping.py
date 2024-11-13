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
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/zoonomia447/hg38.phyloP447way.bw'
phylop_447=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', custom_id=True)
# %%
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/NFkB-p65.bed'
nfkb=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=True)
#%%
mean_phylop=mapper.normalise(alignment=alignment_filtered, mapped_data=phylop_447)
mean_nfkb=mapper.normalise(alignment=alignment_filtered, mapped_data=nfkb)
# %%
#find peaks
import scipy
peaks, _ = scipy.signal.find_peaks(mean_nfkb, width = 6)
#highest_peak_index = peaks[np.argmax(mean_znf808[peaks])]
binding_indices = np.unique(np.where(nfkb[:, peaks] != 0)[0])
nonbinding_indices=list(set(np.arange(nfkb.shape[0])) - set(binding_indices))
#%%
# %%
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
phylop_sorted = np.vstack((phylop_bind,phylop_nonbind))
nfkb_bind = nfkb[binding_indices]
nfkb_nonbind = nfkb[nonbinding_indices]
nfkb_sorted = np.vstack((nfkb_bind,nfkb_nonbind))
te_age_sorted=metadata_age.iloc[np.concatenate((binding_indices,nonbinding_indices))].te_age.fillna(0)
# %%
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')

#%%
anno_label=metadata_age.fillna(0).te_age.sort_values().unique()
# %% 
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [phylop_sorted,nfkb_sorted,], 
    alignment=alignment_filtered, 
    heatmap_color=[custom_cmap.vlag_r_mpl,'Blues'],
    heatmap_mode='spread_horizontal', 
    vlim = [[-0.5,0.5],[0,0.0005],[0,0.0005]], 
    opacity = 0.5, 
    annotation=True, 
    anno_col = ['Blues'], 
    annotation_data=[te_age_sorted],
    anno_cbar_label=[anno_label],
    anno_title=['age'],
    anno_cbar_title=['MYA'], 
    aggregated=True, 
    aggregated_data=[mean_phylop,-np.log10(p_value_greater),np.log10(p_value_less)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[None,None]],
    agg_titles=['mean_phyloP','-log10P p>a','log10P p<a'], 
    colorbar=True,
    colorbar_steps = [0.1,0.0001],  )
#%%