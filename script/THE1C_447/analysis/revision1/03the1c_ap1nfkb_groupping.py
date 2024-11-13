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
from scipy.stats import false_discovery_control
def fdr_control_with_nans(p_values, method='bh'):
    p_values = np.asarray(p_values)
    valid_mask = ~np.isnan(p_values)
    valid_p_values = p_values[valid_mask]
    adj_p_values_valid = false_discovery_control(valid_p_values, method=method)
    adj_p_values = np.full(p_values.shape, np.nan)
    adj_p_values[valid_mask] = adj_p_values_valid
    return adj_p_values
age_ref_table = pd.DataFrame(config.age_ref_table_template)
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
###############IMPORTANT####################
#filter NA
metadata_age=metadata_age[~metadata_age['te_age'].isna()]
noNA_indices=metadata_age.index
metadata_age=metadata_age.reset_index()
alignment_filtered=alignment_filtered[noNA_indices]
#%%
import numpy as np
metadata_age['len'] = metadata_age.end.astype(int) - metadata_age.start.astype(int)
age_subgroups = np.unique(metadata_age['te_age'].sort_values())
age_subgroup = {subgroup: num for num, subgroup in enumerate(age_subgroups)}
age_anno=metadata_age['te_age'].map(age_subgroup)
#%%
coord_file=mapper.extract_metadata_from_alignment(alignment_file)
#%%
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/zoonomia447/hg38.phyloP447way.bw'
phylop_447=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
phylop_447=phylop_447[noNA_indices]
# %%
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/AP-1(bZIP).bed'
ap1=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
ap1=ap1[noNA_indices]
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/NFkB-p65.bed'
nfkb=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
nfkb=nfkb[noNA_indices]
#%%
import importlib
importlib.reload(mapper)
mean_phylop=mapper.normalise(alignment=alignment_filtered, mapped_data=phylop_447)
cov_ap1=mapper.normalise(alignment=alignment_filtered, mapped_data=ap1, method='perc_coverage')
cov_nfkb=mapper.normalise(alignment=alignment_filtered, mapped_data=nfkb, method='perc_coverage')
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [ap1,nfkb], 
    alignment=alignment_filtered,
    show_alignment=True, 
    heatmap_color=['Greens','Blues'],
    heatmap_mode='overlay',
    heatmap_title=['AP-1 and NFkB motifs on THE1C MSA'],
    heatmap_title_fs = 10, 
    anno_ylabel = 'sequences',
    anno_ylabel_fs=10,
    vlim = [[0,7]], 
    opacity = 0.8, 
    hm_transparency_mode='gradient',
    aggregated=True, 
    aggregated_data=[cov_ap1,cov_nfkb], 
    agg_titles=['AP-1','NFkB'], 
    agg_titles_fs=7,
    agg_colset=['green','blue'],
    agg_ylim=[[0,70],[0,70]],
    #agg_titles=['ap1 motif'], 
    agg_ylabel=['perc_coverage', None],
    agg_ylabel_ypos=[0.05,None],
    agg_ylabel_fs=7,
    agg_xlabel='position (bp)',
    annotation=True, 
    anno_col = ['Blues'], 
    annotation_data=[age_anno],
    anno_cbar=True,
    anno_cbar_label=[age_subgroups],
    anno_cbar_title=['TEA-TIME'], 
    colorbar=False,

    )
#%%
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [ap1,nfkb], 
    alignment=alignment_filtered,
    show_alignment=True, 
    heatmap_color=['Greens','Blues'],
    heatmap_mode='overlay',
    heatmap_title=['AP-1 and NFkB motifs on THE1C MSA'],
    heatmap_title_fs = 10, 
    anno_ylabel = 'sequences',
    anno_ylabel_fs=10,
    vlim = [[0,7]], 
    opacity = 0.9,
    hm_transparency_mode = 'gradient', 
    hm_xlabel = 'position (bp)',
    hm_ylabel = 'sequences',
    )
# %%
#find peaks
import scipy
peaks, _ = scipy.signal.find_peaks(cov_ap1, width = 6)
#%%
binding_indices_right = np.unique(np.where(ap1[:, 306] != 0)[0])
nonbinding_indices_right=list(set(np.arange(ap1.shape[0])) - set(binding_indices_right))
#%%
metadata_age['ap1'] = 'No Group'
metadata_age.loc[metadata_age.index.isin(binding_indices_right), 'ap1'] = 'A'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices_right), 'ap1'] = 'B'
# %%
#find peaks
import scipy
peaks, _ = scipy.signal.find_peaks(cov_nfkb, width = 6)
#%%
binding_indices = np.unique(np.where(nfkb[:, peaks] != 0)[0])
nonbinding_indices=list(set(np.arange(nfkb.shape[0])) - set(binding_indices))
#%%
binding_indices_right = np.unique(np.where(nfkb[:, 180] != 0)[0])
nonbinding_indices_right=list(set(np.arange(nfkb.shape[0])) - set(binding_indices_right))

#%%
metadata_age['nfkb'] = 'No Group'
metadata_age.loc[metadata_age.index.isin(binding_indices_right), 'nfkb'] = 'A'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices_right), 'nfkb'] = 'B'

#%%
metadata_sorted=metadata_age.sort_values(['ap1','nfkb','te_age'])
sorted_indices=metadata_sorted.index
subgroups = np.unique(metadata_sorted['ap1'].astype(str))
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno_right=metadata_sorted['ap1'].map(numerical_subgroup)
subgroups = np.unique(metadata_sorted['nfkb'].astype(str))
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno_left=metadata_sorted['nfkb'].map(numerical_subgroup)
#%%
phylop_sorted=phylop_447[sorted_indices]
ap1_sorted=ap1[sorted_indices]
nfkb_sorted=nfkb[sorted_indices]
alignment_sorted=alignment_filtered[sorted_indices]
age_anno_sorted=age_anno[sorted_indices]
#%%

#%%
binding_indices=metadata_sorted[(metadata_sorted['ap1']=='A')|(metadata_sorted['nfkb']=='A')].index
nonbinding_indices=metadata_sorted[(metadata_sorted['ap1']=='B')&(metadata_sorted['nfkb']=='B')].index
only_right = metadata_sorted[(metadata_sorted['ap1']=='A')&(metadata_sorted['nfkb']=='B')].index
only_left = metadata_sorted[(metadata_sorted['ap1']=='B')&(metadata_sorted['nfkb']=='A')].index

intersect_right_left = metadata_sorted[(metadata_sorted['ap1']=='A')&(metadata_sorted['nfkb']=='A')].index

#%%
metadata_age['motif_group'] = 'No Group'
metadata_age.loc[metadata_age.index.isin(intersect_right_left), 'motif_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(only_right), 'motif_group'] = 'B'
metadata_age.loc[metadata_age.index.isin(only_left), 'motif_group'] = 'C'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices), 'motif_group'] = 'D'
subgroups = np.unique(metadata_age['motif_group'].astype(str))
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno=metadata_age['motif_group'].map(numerical_subgroup)
subgroup_anno_sorted=subgroup_anno[sorted_indices]
#%%
#%%
test_set=metadata_sorted.reset_index()
#%%
test_set[(test_set['ap1']=='A')&(test_set['nfkb']=='A')].index
#%%
test_set[(test_set['ap1']=='A')&(test_set['nfkb']=='B')].index
#%%
test_set[(test_set['ap1']=='B')&(test_set['nfkb']=='A')].index
#%%
#%%
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [ap1_sorted,nfkb_sorted,], 
    alignment=alignment_sorted,
    show_alignment=True, 
    heatmap_color=['Greens','Blues',],
    heatmap_mode='overlay', 
    vlim = [[0,7]], 
    opacity = 0.9,
    hm_transparency_mode='gradient', 
    heatmap_title=['AP-1 and NFkB motifs on THE1C, grouped by motifs'],
    heatmap_title_fs = 10, 
    heatmap_xhighlight=[[2624,2624],[5573,5573],[7098,7098]],
    heatmap_xhighlight_col= ['black','black','black'],
    heatmap_xhighlight_alpha=[0.5,0.5,0.5],
    heatmap_xlabel = 'position (bp)',
    heatmap_xlabel_fs=10,
    heatmap_ylabel= 'sequences',
    heatmap_ylabel_fs=10
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
    vlim = [[-1,1]], 
    opacity = 0.8, 
    hm_transparency_mode='gradient',
    hm_interpolation = None,
    annotation=True, 
    anno_col = ['Blues',['purple','red','blue','white']], 
    annotation_data=[age_anno_sorted,subgroup_anno_sorted],
    anno_cbar_label=[age_subgroups,['both','only ap1','only nfkb','no motif/\nothers']],
    anno_cbar_title=['TEA-TIME','motif group'], 
    anno_cbar_even_pos=-0.19,
    colorbar=True,
    colorbar_steps=[0.1],
    heatmap_title=['phyloP of THE1C, grouped by AP-1 and NFkB motifs'],
    heatmap_title_fs = 10, 
    heatmap_xlabel = 'position (bp)',
    heatmap_xlabel_fs=10,
    anno_ylabel = 'sequences',
    anno_ylabel_fs=10,
    heatmap_xhighlight=[[2624,2624],[5573,5573],[7098,7098]],
    heatmap_xhighlight_col= ['black','black','black'],
    heatmap_xhighlight_alpha=[0.8,0.8,0.8],
    )
#%%
metadata_age_filtered=metadata_age
binding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='A')&(metadata_age_filtered['nfkb']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='B')&(metadata_age_filtered['nfkb']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=fdr_control_with_nans(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=fdr_control_with_nans(p_value_less)
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
    aggregated_data=[-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['blue','red'],
    agg_ylim=[[0,20],[-20,0]],
    agg_yhighlight=[[175,190],[304,314]],
    agg_yhighlight_col= ['blue','green'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['conserved\nmotif>no motif','accelerated\nmotif<no motif'], 
    agg_titles_fs=10,
    agg_titles_pos=[1.2,0.5],
    agg_ylabel=['-log10(adjusted p-value)',None],
    agg_ylabel_ypos=[0.05,None],
    agg_ylabel_xpos=0,
    agg_ylabel_fs=10,
    agg_plot_title=['phyloP of THE1C, with both AP-1 and NFkB motifs',None],
    agg_plot_title_fs= 12, 
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    agg_plottext=[None, f'TE w/ motifs: {len(binding_indices)}\nTE w/o motifs: {len(nonbinding_indices)}'],
    agg_plottext_fs=10,
    agg_plottext_pos=[0.99,0.01],
    agg_major_tick=20,
    colorbar=False,
    agg_h=20,
    figsize= [80,40],
    #xlim=[304,313],
    #gg_major_tick=1,
    )
#%%
age_of_interest_list=config.age_canon[0:7]
for idx, age_of_interest in enumerate(age_of_interest_list):
    metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
    binding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='A')&(metadata_age_filtered['nfkb']=='A')].index
    nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='B')&(metadata_age_filtered['nfkb']=='B')].index
    phylop_bind = phylop_447[binding_indices]
    phylop_nonbind = phylop_447[nonbinding_indices]
    alignemnt_bind=alignment_filtered[binding_indices]
    phylop_bind[alignemnt_bind == 0] = np.nan
    alignemnt_nonbind=alignment_filtered[nonbinding_indices]
    phylop_nonbind[alignemnt_nonbind == 0] = np.nan
    alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
    stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
    p_value_greater_adjusted=fdr_control_with_nans(p_value_greater)
    stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
    p_value_less_adjusted=fdr_control_with_nans(p_value_less)
    greater_sig = False
    less_sig = False
    if np.nanmax(-np.log10(p_value_greater_adjusted))>3:
        greater_sig = True  
    if np.nanmax(-np.log10(p_value_less_adjusted))>3:
        less_sig = True
    
    print(f'{age_of_interest}\t{len(binding_indices)}\t{len(nonbinding_indices)}\t{np.nanmin([np.nanmin(p_value_greater_adjusted),np.nanmin(p_value_less_adjusted)]):.2e}')
#%%
age_of_interest = 43.2
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='A')&(metadata_age_filtered['nfkb']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='B')&(metadata_age_filtered['nfkb']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=fdr_control_with_nans(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=fdr_control_with_nans(p_value_less)
#%%
group_name = age_ref_table[age_ref_table['age']==age_of_interest]['representative'].values[0]
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
    aggregated_data=[-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['blue','red'],
    agg_ylim=[[0,20],[-20,0]],
    agg_yhighlight=[[175,190],[304,314]],
    agg_yhighlight_col= ['blue','green'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['conserved\nmotif>no motif','accelerated\nmotif<no motif'], 
    agg_titles_fs=10,
    agg_titles_pos=[1.2,0.5],
    agg_ylabel=['-log10(adjusted p-value)',None],
    agg_ylabel_ypos=[0.05,None],
    agg_ylabel_xpos=0,
    agg_ylabel_fs=10,
    #agg_plot_title=['phyloP of THE1C, grouped by AP-1 motif',None],
    agg_plot_title_fs= 12, 
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    agg_plottext=[None,f'{group_name} ({age_of_interest} MYA)\nTE w/ motifs: {len(binding_indices)}\nTE w/o motifs: {len(nonbinding_indices)}'],
    agg_plottext_fs=10,
    agg_plottext_pos=[0.99,0.01],
    agg_major_tick=20,
    colorbar=False,
    agg_h=20,
    figsize= [80,40],
    #xlim=[304,313],
    #gg_major_tick=1,
    )
#%%
metadata_age_filtered=metadata_age
binding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='A')&(metadata_age_filtered['nfkb']=='B')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='B')&(metadata_age_filtered['nfkb']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=fdr_control_with_nans(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=fdr_control_with_nans(p_value_less)
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
    aggregated_data=[-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['blue','red'],
    agg_ylim=[[0,15],[-15,0]],
    agg_yhighlight=[[175,190],[304,314]],
    agg_yhighlight_col= ['blue','green'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['conserved\nmotif>no motif','accelerated\nmotif<no motif'], 
    agg_titles_fs=10,
    agg_titles_pos=[1.2,0.5],
    agg_ylabel=['-log10(adjusted p-value)',None],
    agg_ylabel_ypos=[0.05,None],
    agg_ylabel_xpos=0,
    agg_ylabel_fs=10,
    agg_plot_title=['phyloP of THE1C, with only AP-1 motif',None],
    agg_plot_title_fs= 12, 
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    agg_plottext=[None, f'TE w/ motifs: {len(binding_indices)}\nTE w/o motifs: {len(nonbinding_indices)}'],
    agg_plottext_fs=10,
    agg_plottext_pos=[0.99,0.01],
    agg_major_tick=20,
    colorbar=False,
    agg_h=20,
    figsize= [80,40],
    #xlim=[304,313],
    #gg_major_tick=1,
    )
#%%
age_of_interest_list=config.age_canon[0:7]
for idx, age_of_interest in enumerate(age_of_interest_list):
    metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
    binding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='A')&(metadata_age_filtered['nfkb']=='B')].index
    nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='B')&(metadata_age_filtered['nfkb']=='B')].index
    phylop_bind = phylop_447[binding_indices]
    phylop_nonbind = phylop_447[nonbinding_indices]
    alignemnt_bind=alignment_filtered[binding_indices]
    phylop_bind[alignemnt_bind == 0] = np.nan
    alignemnt_nonbind=alignment_filtered[nonbinding_indices]
    phylop_nonbind[alignemnt_nonbind == 0] = np.nan
    alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
    stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
    p_value_greater_adjusted=fdr_control_with_nans(p_value_greater)
    stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
    p_value_less_adjusted=fdr_control_with_nans(p_value_less)
    greater_sig = False
    less_sig = False
    if np.nanmax(-np.log10(p_value_greater_adjusted))>3:
        greater_sig = True  
    if np.nanmax(-np.log10(p_value_less_adjusted))>3:
        less_sig = True
    
    print(f'{age_of_interest}\t{len(binding_indices)}\t{len(nonbinding_indices)}\t{np.nanmin([np.nanmin(p_value_greater_adjusted),np.nanmin(p_value_less_adjusted)]):.2e}')
#%%
age_of_interest = 43.2
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='A')&(metadata_age_filtered['nfkb']=='B')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='B')&(metadata_age_filtered['nfkb']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=fdr_control_with_nans(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=fdr_control_with_nans(p_value_less)
#%%
group_name = age_ref_table[age_ref_table['age']==age_of_interest]['representative'].values[0]
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
    aggregated_data=[-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['blue','red'],
    agg_ylim=[[0,15],[-15,0]],
    agg_yhighlight=[[175,190],[304,314]],
    agg_yhighlight_col= ['blue','green'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['conserved\nmotif>no motif','accelerated\nmotif<no motif'], 
    agg_titles_fs=10,
    agg_titles_pos=[1.2,0.5],
    agg_ylabel=['-log10(adjusted p-value)',None],
    agg_ylabel_ypos=[0.05,None],
    agg_ylabel_xpos=0,
    agg_ylabel_fs=10,
    #agg_plot_title=['phyloP of THE1C, grouped by AP-1 motif',None],
    agg_plot_title_fs= 12, 
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    agg_plottext=[None,f'{group_name} ({age_of_interest} MYA)\nTE w/ motifs: {len(binding_indices)}\nTE w/o motifs: {len(nonbinding_indices)}'],
    agg_plottext_fs=10,
    agg_plottext_pos=[0.99,0.01],
    agg_major_tick=20,
    colorbar=False,
    agg_h=20,
    figsize= [80,40],
    #xlim=[304,313],
    #gg_major_tick=1,
    )
#%%
metadata_age_filtered=metadata_age
binding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='B')&(metadata_age_filtered['nfkb']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='B')&(metadata_age_filtered['nfkb']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=fdr_control_with_nans(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=fdr_control_with_nans(p_value_less)
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
    aggregated_data=[-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['blue','red'],
    agg_ylim=[[0,10],[-10,0]],
    agg_yhighlight=[[175,190],[304,314]],
    agg_yhighlight_col= ['blue','green'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['conserved\nmotif>no motif','accelerated\nmotif<no motif'], 
    agg_titles_fs=10,
    agg_titles_pos=[1.2,0.5],
    agg_ylabel=['-log10(adjusted p-value)',None],
    agg_ylabel_ypos=[0.05,None],
    agg_ylabel_xpos=0,
    agg_ylabel_fs=10,
    agg_plot_title=['phyloP of THE1C, with only NFkB motif',None],
    agg_plot_title_fs= 12, 
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    agg_plottext=[None, f'TE w/ motifs: {len(binding_indices)}\nTE w/o motifs: {len(nonbinding_indices)}'],
    agg_plottext_fs=10,
    agg_plottext_pos=[0.99,0.01],
    agg_major_tick=20,
    colorbar=False,
    agg_h=20,
    figsize= [80,40],
    #xlim=[304,313],
    #gg_major_tick=1,
    )
#%%
age_of_interest_list=config.age_canon[0:7]
for idx, age_of_interest in enumerate(age_of_interest_list):
    metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
    binding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='B')&(metadata_age_filtered['nfkb']=='A')].index
    nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='B')&(metadata_age_filtered['nfkb']=='B')].index
    phylop_bind = phylop_447[binding_indices]
    phylop_nonbind = phylop_447[nonbinding_indices]
    alignemnt_bind=alignment_filtered[binding_indices]
    phylop_bind[alignemnt_bind == 0] = np.nan
    alignemnt_nonbind=alignment_filtered[nonbinding_indices]
    phylop_nonbind[alignemnt_nonbind == 0] = np.nan
    alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
    stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
    p_value_greater_adjusted=fdr_control_with_nans(p_value_greater)
    stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
    p_value_less_adjusted=fdr_control_with_nans(p_value_less)
    greater_sig = False
    less_sig = False
    if np.nanmax(-np.log10(p_value_greater_adjusted))>3:
        greater_sig = True  
    if np.nanmax(-np.log10(p_value_less_adjusted))>3:
        less_sig = True
    
    print(f'{age_of_interest}\t{len(binding_indices)}\t{len(nonbinding_indices)}\t{np.nanmin([np.nanmin(p_value_greater_adjusted),np.nanmin(p_value_less_adjusted)]):.2e}')
#%%
age_of_interest = 43.2
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='B')&(metadata_age_filtered['nfkb']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='B')&(metadata_age_filtered['nfkb']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=fdr_control_with_nans(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=fdr_control_with_nans(p_value_less)
#%%
group_name = age_ref_table[age_ref_table['age']==age_of_interest]['representative'].values[0]
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
    aggregated_data=[-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['blue','red'],
    agg_ylim=[[0,10],[-10,0]],
    agg_yhighlight=[[175,190],[304,314]],
    agg_yhighlight_col= ['blue','green'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['conserved\nmotif>no motif','accelerated\nmotif<no motif'], 
    agg_titles_fs=10,
    agg_titles_pos=[1.2,0.5],
    agg_ylabel=['-log10(adjusted p-value)',None],
    agg_ylabel_ypos=[0.05,None],
    agg_ylabel_xpos=0,
    agg_ylabel_fs=10,
    #agg_plot_title=['phyloP of THE1C, grouped by AP-1 motif',None],
    agg_plot_title_fs= 12, 
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    agg_plottext=[None,f'{group_name} ({age_of_interest} MYA)\nTE w/ motifs: {len(binding_indices)}\nTE w/o motifs: {len(nonbinding_indices)}'],
    agg_plottext_fs=10,
    agg_plottext_pos=[0.99,0.01],
    agg_major_tick=20,
    colorbar=False,
    agg_h=20,
    figsize= [80,40],
    #xlim=[304,313],
    #gg_major_tick=1,
    )

#%%
#%%
age_of_interest = 20.19
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='B')&(metadata_age_filtered['nfkb']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['ap1']=='B')&(metadata_age_filtered['nfkb']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=fdr_control_with_nans(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=fdr_control_with_nans(p_value_less)
#%%
group_name = age_ref_table[age_ref_table['age']==age_of_interest]['representative'].values[0]
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
    aggregated_data=[-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['blue','red'],
    agg_ylim=[[0,10],[-10,0]],
    agg_yhighlight=[[175,190],[304,314]],
    agg_yhighlight_col= ['blue','green'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['conserved\nmotif>no motif','accelerated\nmotif<no motif'], 
    agg_titles_fs=10,
    agg_titles_pos=[1.2,0.5],
    agg_ylabel=['-log10(adjusted p-value)',None],
    agg_ylabel_ypos=[0.05,None],
    agg_ylabel_xpos=0,
    agg_ylabel_fs=10,
    #agg_plot_title=['phyloP of THE1C, grouped by AP-1 motif',None],
    agg_plot_title_fs= 12, 
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    agg_plottext=[None,f'{group_name} ({age_of_interest} MYA)\nTE w/ motifs: {len(binding_indices)}\nTE w/o motifs: {len(nonbinding_indices)}'],
    agg_plottext_fs=10,
    agg_plottext_pos=[0.99,0.01],
    agg_major_tick=20,
    colorbar=False,
    agg_h=20,
    figsize= [80,40],
    #xlim=[304,313],
    #gg_major_tick=1,
    )
# %%
