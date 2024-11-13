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
subfamily='MER11'
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
#%%
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
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
import numpy as np
metadata_age['len'] = metadata_age.end.astype(int) - metadata_age.start.astype(int)
metadata_age['subgroup'] = metadata_age['group'].str.split('_').str[0]
subgroups = np.unique(metadata_age['subgroup'].astype(str))
#%%
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno=metadata_age['subgroup'].map(numerical_subgroup)
# %%
coord_file = mapper.extract_metadata_from_alignment(alignment_file)
#%%

vcf_dir = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/vcf-gnomad'

vcf=mapper.map_and_overlay(alignment_file, coord_file, vcf_dir, data_format='vcf', vcf_format='gnomad', custom_id=True)
# %%
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/kzfp_peak_bed/znf808_lifted.bed'
znf808=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap = False)
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/kzfp_peak_bed/znf525_lifted.bed'
znf525=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap = False)

#%%
mean_vcf=mapper.normalise(alignment=alignment_filtered, mapped_data=vcf)
mean_znf808=mapper.normalise(alignment=alignment_filtered, mapped_data=znf808)
mean_znf525=mapper.normalise(alignment=alignment_filtered, mapped_data=znf525)
# %%
#find peaks
import scipy
peaks, _ = scipy.signal.find_peaks(mean_znf808, width = 8)
#highest_peak_index = peaks[np.argmax(mean_znf808[peaks])]
binding_indices = np.unique(np.where(znf808[:, peaks] != 0)[0])
nonbinding_indices=list(set(np.arange(znf808.shape[0])) - set(binding_indices))
# %%
vcf_bind = vcf[binding_indices]
vcf_nonbind = vcf[nonbinding_indices]
vcf_sorted = np.vstack((vcf_bind,vcf_nonbind))
znf808_bind = znf808[binding_indices]
znf808_nonbind = znf808[nonbinding_indices]
znf808_sorted = np.vstack((znf808_bind,znf808_nonbind))
znf525_bind = znf525[binding_indices]
znf525_nonbind = znf525[nonbinding_indices]
znf525_sorted = np.vstack((znf525_bind,znf525_nonbind))
te_age_sorted=metadata_age.iloc[np.concatenate((binding_indices,nonbinding_indices))].te_age.fillna(0)
subgroup_anno_sorted=subgroup_anno[np.concatenate((binding_indices,nonbinding_indices))]
# %%
alignemnt_bind=alignment_filtered[binding_indices]
vcf_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
vcf_nonbind[alignemnt_nonbind == 0] = np.nan
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(vcf_bind,vcf_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(vcf_bind,vcf_nonbind, axis =0,nan_policy='omit', alternative = 'less')

# %% experimental - spread aggregated  and heatmap
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [vcf_sorted,znf808_sorted,znf525_sorted], alignment=alignment_filtered, 
    heatmap_color=['Blues','Greens','Oranges'],
    heatmap_title=['MER11A AAF sorted by ZNF808 binding','ChIP-exo ZNF808','ChIP-exo ZNF525'],  
    vlim = [[0,0.0001],[0,0.0005],[0,0.0005]], 
    opacity = 0.5,
    heatmap_mode='spread_horizontal',  
    annotation=True, 
    annotation_data=[te_age_sorted, subgroup_anno_sorted],
    anno_col = ['Blues',['red', 'yellow', 'blue']], 
    anno_title=['age','subfamily'], 
    anno_cbar_title=['age (MYA)','subfamily'], 
    anno_cbar_label=[te_age_sorted.sort_values().unique(),['MER11A','MER11B','MER11C']],
    aggregated=True, 
    aggregated_data=[mean_vcf,-np.log10(p_value_greater),np.log10(p_value_less)], 
    agg_colset=['grey','blue','red'], 
    agg_ylim = [[None,None], [None,None], [None,None]], 
    agg_ylabel=['AAF','-log10P','log10P'], 
    agg_major_tick=100,
    agg_titles=['mean AAF','ZNF808 p>a', 'ZNF808 p<a'],
    colorbar=True, 
    colorbar_steps=[0.00001,10,10], 
    )
#%% filter by age
#29.44 20.19 15.76 9.06 6.7 
metadata_age_subset=metadata_age[metadata_age.te_age==6.7]
alignment_filtered_subset = alignment_filtered[metadata_age_subset.index]
subgroups_subset = np.unique(metadata_age_subset['subgroup'].astype(str))
numerical_subgroup_subset = {subgroups_subset: num for num, subgroups_subset in enumerate(subgroups_subset)}
subgroup_anno_subset=metadata_age_subset['subgroup'].map(numerical_subgroup_subset)
coord_subset=coord_file.iloc[metadata_age_subset.index]
znf808_subset=znf808[metadata_age_subset.index]
vcf_subset=vcf[metadata_age_subset.index]
# %%
mean_vcf=mapper.normalise(alignment=alignment_filtered_subset, mapped_data=vcf_subset)
mean_znf808=mapper.normalise(alignment=alignment_filtered_subset, mapped_data=znf808_subset)

#find peaks
import scipy
peaks, _ = scipy.signal.find_peaks(mean_znf808, width = 8)
#highest_peak_index = peaks[np.argmax(mean_znf808[peaks])]
binding_indices = np.unique(np.where(znf808_subset[:, peaks] != 0)[0])
nonbinding_indices=list(set(np.arange(znf808_subset.shape[0])) - set(binding_indices))
#%%
te_age_sorted_subset=metadata_age_subset.reset_index().iloc[np.concatenate((binding_indices,nonbinding_indices))].te_age.fillna(0)
#%%
subgroup_anno_sorted_subset=subgroup_anno.iloc[np.concatenate((binding_indices,nonbinding_indices))]
vcf_bind_subset = vcf_subset[binding_indices]
vcf_nonbind_subset = vcf_subset[nonbinding_indices]
vcf_sorted = np.vstack((vcf_bind_subset,vcf_nonbind_subset))
znf808_bind_subset = znf808_subset[binding_indices]
znf808_nonbind_subset = znf808_subset[nonbinding_indices]
znf808_sorted = np.vstack((znf808_bind_subset,znf808_nonbind_subset))
# %%
alignemnt_bind=alignment_filtered_subset[binding_indices]
vcf_bind_subset[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered_subset[nonbinding_indices]
vcf_nonbind_subset[alignemnt_nonbind == 0] = np.nan
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(vcf_bind_subset,vcf_nonbind_subset, axis =0,nan_policy='omit', alternative = 'greater')
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(vcf_bind_subset,vcf_nonbind_subset, axis =0,nan_policy='omit', alternative = 'less')
#%%
# %% experimental - spread aggregated  and heatmap
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [vcf_sorted,znf808_sorted], 
    alignment=alignment_filtered_subset, 
    heatmap_mode='spread_horizontal',
    heatmap_color=['Blues','Greens'],
    heatmap_title=['MER11A AAF sorted by ZNF808 binding','ChIP-exo ZNF808'], 
    vlim = [[0,0.0001],[0,0.0005]], 
    opacity = 0.5, 
    annotation=True, 
    annotation_data=[ subgroup_anno_subset],
    anno_col = [['red', 'yellow', 'blue']],
    anno_cbar_title=['subfamily'], 
    anno_cbar_label=[['MER11A','MER11B','MER11C']],
    anno_title=['subfamily'],
    aggregated=True, 
    aggregated_data=[mean_vcf,-np.log10(p_value_greater),np.log10(p_value_less)], 
    agg_colset=['grey','blue','red'],
    agg_titles=['mean AAF','ZNF808 p>a','ZNF808 p<a'], 
    agg_ylabel=['AAF','-log10P','log10P'],
    agg_ylim = [[None,None], [None,None],[None,None]],  
    agg_major_tick=100, 
    colorbar=True, 
    colorbar_steps=[0.1,10,10])
# %%
# debug session 
from ma_mapper import extract_bed
metadata_bed=extract_bed.metadata_to_bed(coord_file)
import pybedtools
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/kzfp_peak_bed/znf808_lifted.bed'
bed_file = pybedtools.BedTool(bed_file)
metadata_bed=pybedtools.BedTool.from_dataframe(metadata_bed.iloc[:,0:6])
intersect_bed=metadata_bed.intersect(bed_file,loj=True, wa = True, wb =True, s=False)
intersect_df=intersect_bed.to_dataframe()
intersect_df[intersect_df['thickEnd']!=-1]
# %%
