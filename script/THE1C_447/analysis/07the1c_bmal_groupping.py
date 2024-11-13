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
subfamily = 'THE1C'
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file,col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
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
phylop_447=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
# %%
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/BMAL1.bed'
bmal=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
#%%
import importlib
importlib.reload(mapper)
mean_phylop=mapper.normalise(alignment=alignment_filtered, mapped_data=phylop_447)
cov_bmal=mapper.normalise(alignment=alignment_filtered, mapped_data=bmal, method='perc_coverage')
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [bmal,], 
    alignment=alignment_filtered,
    show_alignment=True, 
    heatmap_color=['Reds'],
    heatmap_mode='overlay', 
    vlim = [[0,7]], 
    opacity = 0.5, 
    aggregated=True, 
    aggregated_data=[cov_bmal], 
    agg_colset=['red',],
    agg_ylim=[[0,70]],
    agg_titles=['BMAL1 motif'], 
    agg_ylabel=['perc_coverage'],
    colorbar=False,
    )
# %%
#find peaks
import scipy
peaks, _ = scipy.signal.find_peaks(cov_bmal, width = 6)
#%%
#65,227/240,294
#highest_peak_index = peaks[np.argmax(mean_znf808[peaks])]
binding_indices = np.unique(np.where(bmal[:, peaks] != 0)[0])
nonbinding_indices=list(set(np.arange(bmal.shape[0])) - set(binding_indices))
#%%
binding_indices_65 = np.unique(np.where(bmal[:, 65] != 0)[0])
nonbinding_indices_65=list(set(np.arange(bmal.shape[0])) - set(binding_indices_65))
binding_indices_227 = np.unique(np.where(bmal[:, 227] != 0)[0])
nonbinding_indices_227=list(set(np.arange(bmal.shape[0])) - set(binding_indices_227))
binding_indices_294 = np.unique(np.where(bmal[:, 294] != 0)[0])
nonbinding_indices_294=list(set(np.arange(bmal.shape[0])) - set(binding_indices_294))
#%%
metadata_age['65_group'] = 'No Group'
metadata_age.loc[metadata_age.index.isin(binding_indices_65), '65_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices_65), '65_group'] = 'B'
metadata_age['227_group'] = 'No Group'
metadata_age.loc[metadata_age.index.isin(binding_indices_227), '227_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices_227), '227_group'] = 'B'
metadata_age['294_group'] = 'No Group'
metadata_age.loc[metadata_age.index.isin(binding_indices_294), '294_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices_294), '294_group'] = 'B'

#%%
metadata_sorted=metadata_age.sort_values(['65_group','227_group','294_group'])
final_sort_indices=metadata_sorted.index
subgroups = np.unique(metadata_sorted['65_group'].astype(str))
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno_65=metadata_sorted['65_group'].map(numerical_subgroup)
subgroups = np.unique(metadata_sorted['227_group'].astype(str))
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno_227=metadata_sorted['227_group'].map(numerical_subgroup)
subgroups = np.unique(metadata_sorted['294_group'].astype(str))
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno_294=metadata_sorted['294_group'].map(numerical_subgroup)

#%%

#cov_bmal[61:69]
#array([31.10003477, 60.56826057, 60.95469818, 61.43138847, 61.60027952,
#       61.5608281 , 61.18610823, 30.90846207])
#cov_bmal[224:232]
#array([23.36633663, 46.68298939, 46.72494172, 47.47749882, 46.95449242,
#       46.96545285, 46.89509999, 24.74657881])
#cov_bmal[291:299]
#array([ 6.24419684, 12.77906977, 12.82852842, 12.89194545, 12.80713955,
#       12.78813854, 12.74441341,  6.54622243])
#%%
phylop_sorted=phylop_447[final_sort_indices]
bmal_sorted=bmal[final_sort_indices]
alignment_sorted=alignment_filtered[final_sort_indices]
#%%

#%%
binding_indices=metadata_sorted[(metadata_sorted['65_group']=='A')|(metadata_sorted['227_group']=='A')|(metadata_sorted['294_group']=='A')].index
nonbinding_indices=metadata_sorted[(metadata_sorted['65_group']=='B')&(metadata_sorted['227_group']=='B')&(metadata_sorted['294_group']=='B')].index
only_65 = metadata_sorted[(metadata_sorted['65_group']=='A')&(metadata_sorted['227_group']=='B')&(metadata_sorted['294_group']=='B')].index
only_227 = metadata_sorted[(metadata_sorted['65_group']=='B')&(metadata_sorted['227_group']=='A')&(metadata_sorted['294_group']=='B')].index
only_294 = metadata_sorted[(metadata_sorted['65_group']=='B')&(metadata_sorted['227_group']=='B')&(metadata_sorted['294_group']=='A')].index
intersect_65_227 = metadata_sorted[(metadata_sorted['65_group']=='A')&(metadata_sorted['227_group']=='A')&(metadata_sorted['294_group']=='B')].index
intersect_65_294 = metadata_sorted[(metadata_sorted['65_group']=='A')&(metadata_sorted['227_group']=='B')&(metadata_sorted['294_group']=='A')].index
intersect_227_294 = metadata_sorted[(metadata_sorted['65_group']=='B')&(metadata_sorted['227_group']=='A')&(metadata_sorted['294_group']=='A')].index
intersect_65_227_294 = metadata_sorted[(metadata_sorted['65_group']=='A')&(metadata_sorted['227_group']=='A')&(metadata_sorted['294_group']=='A')].index
#%%
metadata_age['motif_group'] = 'No Group'
metadata_age.loc[metadata_age.index.isin(intersect_65_227_294), 'motif_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(intersect_65_227), 'motif_group'] = 'B'
metadata_age.loc[metadata_age.index.isin(intersect_65_294), 'motif_group'] = 'C'
metadata_age.loc[metadata_age.index.isin(only_65), 'motif_group'] = 'D'
metadata_age.loc[metadata_age.index.isin(intersect_227_294), 'motif_group'] = 'E'
metadata_age.loc[metadata_age.index.isin(only_227), 'motif_group'] = 'F'
metadata_age.loc[metadata_age.index.isin(only_294), 'motif_group'] = 'G'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices), 'motif_group'] = 'H'
subgroups = np.unique(metadata_age['motif_group'].astype(str))
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno=metadata_age['motif_group'].map(numerical_subgroup)
subgroup_anno_sorted=subgroup_anno.reindex(final_sort_indices)
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [bmal_sorted], 
    alignment=alignment_sorted,
    show_alignment=True, 
    heatmap_color=['Reds'],
    heatmap_mode='overlay', 
    vlim = [[0,7]], 
    opacity = 0.5, 

    )
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [phylop_sorted,], 
    alignment=alignment_sorted,
    #show_alignment=True, 
    heatmap_color=[custom_cmap.vlag_r_mpl],
    heatmap_mode='overlay', 
    vlim = [[-0.5,0.5]], 
    opacity = 0.5, 
    annotation=True, 
    anno_col = [['black','purple','orange','red','green','blue','yellow','white']], 
    annotation_data=[subgroup_anno_sorted],
    anno_cbar_label=[['65+227+294','65+227','65+294','65 only','227+294','227 only','294 only','no motif/others']],
    anno_cbar_title=['motif group'], 
    colorbar=True,
    colorbar_steps=[0.1]
    )
#%%

phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(p_value_less)
#%%
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    heatmap=False,
    show_alignment=False,
    alignment=alignment_sorted,
    #logos=True, 
    aggregated=True, 
    aggregated_data=[mean_phylop,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,5],[-5,0]],
    agg_yhighlight=[[61,69],[224,232],[291,299]],
    agg_yhighlight_col= ['red','red','red'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    #xlim=[304,313],
    #gg_major_tick=1,
    )
#%%
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    heatmap=False,
    show_alignment=False,
    alignment=alignment_sorted,
    logos=True, 
    aggregated=True, 
    aggregated_data=[mean_phylop,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,5],[-5,0]],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    xlim=[61,69],
    agg_major_tick=1,
    )

# %%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    heatmap=False,
    show_alignment=False,
    alignment=alignment_sorted,
    logos=True, 
    aggregated=True, 
    aggregated_data=[mean_phylop,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,5],[-5,0]],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    xlim=[224,232],
    agg_major_tick=1,
    )
# %%
phylop_bind = phylop_447[intersect_65_227_294]
phylop_nonbind = phylop_447[nonbinding_indices]

alignemnt_bind=alignment_filtered[intersect_65_227_294]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(p_value_less)
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    heatmap=False,
    show_alignment=False, 
    aggregated=True, 
    aggregated_data=[mean_phylop,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[61,69],[224,232],[291,299]],
    agg_yhighlight_col= ['red','red','red'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    agg_h=5
    )
# %%
phylop_bind = phylop_447[intersect_65_227]
phylop_nonbind = phylop_447[nonbinding_indices]

alignemnt_bind=alignment_filtered[intersect_65_227]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(p_value_less)
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    heatmap=False,
    show_alignment=False, 
    aggregated=True, 
    aggregated_data=[mean_phylop,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,10],[-15,0]],
    agg_yhighlight=[[61,69],[224,232],[291,299]],
    agg_yhighlight_col= ['red','red','red'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    agg_h=5
    )

# %%
phylop_bind = phylop_447[intersect_65_294]
phylop_nonbind = phylop_447[nonbinding_indices]

alignemnt_bind=alignment_filtered[intersect_65_294]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(p_value_less)
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    heatmap=False,
    show_alignment=False, 
    aggregated=True, 
    aggregated_data=[mean_phylop,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[61,69],[224,232],[291,299]],
    agg_yhighlight_col= ['red','red','red'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    agg_h=5
    )
# %%
phylop_bind = phylop_447[intersect_227_294]
phylop_nonbind = phylop_447[nonbinding_indices]

alignemnt_bind=alignment_filtered[intersect_227_294]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(p_value_less)
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    heatmap=False,
    show_alignment=False, 
    aggregated=True, 
    aggregated_data=[mean_phylop,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[61,69],[224,232],[291,299]],
    agg_yhighlight_col= ['red','red','red'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    agg_h=5
    )
#%%
phylop_bind = phylop_447[only_65]
phylop_nonbind = phylop_447[nonbinding_indices]

alignemnt_bind=alignment_filtered[only_65]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(p_value_less)
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    heatmap=False,
    show_alignment=False, 
    aggregated=True, 
    aggregated_data=[mean_phylop,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,5],[-5,0]],
    agg_yhighlight=[[61,69],[224,232],[291,299]],
    agg_yhighlight_col= ['red','red','red'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    agg_h=5
    )
#%%
phylop_bind = phylop_447[only_227]
phylop_nonbind = phylop_447[nonbinding_indices]

alignemnt_bind=alignment_filtered[only_227]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(p_value_less)
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    heatmap=False,
    show_alignment=False, 
    aggregated=True, 
    aggregated_data=[mean_phylop,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,10],[-15,0]],
    agg_yhighlight=[[61,69],[224,232],[291,299]],
    agg_yhighlight_col= ['red','red','red'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    agg_h=5
    )
# %%
phylop_bind = phylop_447[only_294]
phylop_nonbind = phylop_447[nonbinding_indices]

alignemnt_bind=alignment_filtered[only_294]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(p_value_less)
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    heatmap=False,
    show_alignment=False, 
    aggregated=True, 
    aggregated_data=[mean_phylop,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[61,69],[224,232],[291,299]],
    agg_yhighlight_col= ['red','red','red'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    agg_h=5
    )
# %%
