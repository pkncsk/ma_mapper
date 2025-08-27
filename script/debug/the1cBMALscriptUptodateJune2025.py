#%% LOAD PACKAGE
from ma_mapper import mapper
from ma_mapper import plots
from ma_mapper import custom_cmap
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
def fdr_control_with_nans(p_values, method='bh'):
    p_values = np.asarray(p_values)
    valid_mask = ~np.isnan(p_values)
    valid_p_values = p_values[valid_mask]
    adj_p_values_valid = false_discovery_control(valid_p_values, method=method)
    adj_p_values = np.full(p_values.shape, np.nan)
    adj_p_values[valid_mask] = adj_p_values_valid
    return adj_p_values
#%% INITIAL PARAMETER
subfamily='THE1C'
alignment_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/alignment/THE1C.fasta.aligned'
genomewide_data_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/annotation/homer_known_motif_hg38/BMAL1(bHLH).bed'
phyloP_data_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/UCSC_phyloP_track/hg38.phyloP447way.bw'
te_age_folder = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime/age_lenient/'
age_table = f'{te_age_folder}/{subfamily}.txt'
internal_id_folder = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime/internal_id'
internal_id_tbl = f'{internal_id_folder}/{subfamily}.internal_id.txt'
#%% SIMPLE WORKFLOW
filtered_alignment_matrix, alignment_coordinate  = mapper.parse_and_filter(alignment_file=alignment_filepath,col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
motif_matrix=mapper.map_and_overlay(alignment_filepath, genomewide_data_filepath,data_format='bed', pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
phyloP_matrix=mapper.map_and_overlay(alignment_filepath, phyloP_data_filepath,data_format='bigwig', pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
age_df = pd.read_csv(age_table, sep='\t')
internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
internal_id_sort = internal_id_df.sort_values('rmsk_index')
te_age_internal_id=internal_id_sort.merge(age_df, on='internal_id', how='left')
#%%
age_default_id = pd.DataFrame()
age_default_id['internal_id'] = subfamily + '_' + te_age_internal_id.index.astype(str)
age_default_id['te_age'] = te_age_internal_id['te_age']
#%%
coord_age = mapper.match_age_to_id_coordinate(alignment_coordinate, age_table=age_default_id)
###############IMPORTANT####################
#filter NA
coord_age=coord_age[~coord_age['te_age'].isna()]
noNA_indices=coord_age.index
coord_age=coord_age.reset_index()
filtered_alignment=filtered_alignment_matrix[noNA_indices]

coord_age['len'] = coord_age.end.astype(int) - coord_age.start.astype(int)
age_subgroups = np.unique(coord_age['te_age'].sort_values())
age_subgroup = {subgroup: num for num, subgroup in enumerate(age_subgroups)}
age_anno=coord_age['te_age'].map(age_subgroup)

phyloP_matrix=phyloP_matrix[noNA_indices]
motif_matrix=motif_matrix[noNA_indices]
# %%
mean_phylop=mapper.normalise(alignment_matrix=filtered_alignment, data_matrix=phyloP_matrix)
cov_motif_matrix=mapper.normalise(alignment_matrix=filtered_alignment, data_matrix=motif_matrix, method='perc_coverage')
#%%
from ma_mapper import plots
from ma_mapper import mapper
plots.plot(
    data = [motif_matrix,], 
    alignment=filtered_alignment,
    show_alignment=True, 
    heatmap_color=['Oranges'],
    heatmap_mode='overlay',
    heatmap_title=['CLOCK motifs on THE1C MSA'],
    heatmap_title_fs = 10, 
    anno_ylabel = 'sequences',
    anno_ylabel_fs=10,
    vlim = [[0,1]], 
    opacity = 0.9, 
    hm_transparency_mode= 'gradient',
    aggregated=True, 
    aggregated_data=[cov_motif_matrix], 
    agg_colset=['orange',],
    agg_ylim=[[0,25]],
    #agg_ylabel_right=['ap1 motif'], 
    agg_ylabel=['perc_coverage'],
    agg_ylabel_fs=7,
    agg_xlabel='position (bp)',
    annotation=False, 
    anno_col = ['Blues'], 
    annotation_data=[age_anno.values],
    anno_cbar=True,
    anno_cbar_label=[age_subgroups],
    anno_cbar_title=['TEA-TIME'], 
    colorbar=False,

    )
# %%
import scipy
peaks, _ = scipy.signal.find_peaks(cov_motif_matrix, width = 6)
#%%
#right,left
#highest_peak_index = peaks[np.argmax(mean_znf808[peaks])]
binding_indices = np.unique(np.where(motif_matrix[:, peaks] != 0)[0])
nonbinding_indices=list(set(np.arange(motif_matrix.shape[0])) - set(binding_indices))
#%%
binding_indices_65 = np.unique(np.where(motif_matrix[:, 65] != 0)[0])
nonbinding_indices_65=list(set(np.arange(motif_matrix.shape[0])) - set(binding_indices_65))
binding_indices_227 = np.unique(np.where(motif_matrix[:, 227] != 0)[0])
nonbinding_indices_227=list(set(np.arange(motif_matrix.shape[0])) - set(binding_indices_227))
binding_indices_294 = np.unique(np.where(motif_matrix[:, 294] != 0)[0])
nonbinding_indices_294=list(set(np.arange(motif_matrix.shape[0])) - set(binding_indices_294))
#%%
coord_age['65_group'] = 'No Group'
coord_age.loc[coord_age.index.isin(binding_indices_65), '65_group'] = 'A'
coord_age.loc[coord_age.index.isin(nonbinding_indices_65), '65_group'] = 'B'
coord_age['227_group'] = 'No Group'
coord_age.loc[coord_age.index.isin(binding_indices_227), '227_group'] = 'A'
coord_age.loc[coord_age.index.isin(nonbinding_indices_227), '227_group'] = 'B'
coord_age['294_group'] = 'No Group'
coord_age.loc[coord_age.index.isin(binding_indices_294), '294_group'] = 'A'
coord_age.loc[coord_age.index.isin(nonbinding_indices_294), '294_group'] = 'B'

#%%
coord_sorted=coord_age.sort_values(['65_group','227_group','294_group'])
final_sort_indices=coord_sorted.index
subgroups = np.unique(coord_sorted['65_group'].astype(str))
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno_65=coord_sorted['65_group'].map(numerical_subgroup)
subgroups = np.unique(coord_sorted['227_group'].astype(str))
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno_227=coord_sorted['227_group'].map(numerical_subgroup)
subgroups = np.unique(coord_sorted['294_group'].astype(str))
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno_294=coord_sorted['294_group'].map(numerical_subgroup)
#%%
phylop_sorted=phyloP_matrix[final_sort_indices]
motif_matrix_sorted=motif_matrix[final_sort_indices]
alignment_sorted=filtered_alignment[final_sort_indices]
age_anno_sorted=age_anno[final_sort_indices]
#%%
binding_indices=coord_sorted[(coord_sorted['65_group']=='A')|(coord_sorted['227_group']=='A')|(coord_sorted['294_group']=='A')].index
nonbinding_indices=coord_sorted[(coord_sorted['65_group']=='B')&(coord_sorted['227_group']=='B')&(coord_sorted['294_group']=='B')].index
only_65 = coord_sorted[(coord_sorted['65_group']=='A')&(coord_sorted['227_group']=='B')&(coord_sorted['294_group']=='B')].index
only_227 = coord_sorted[(coord_sorted['65_group']=='B')&(coord_sorted['227_group']=='A')&(coord_sorted['294_group']=='B')].index
only_294 = coord_sorted[(coord_sorted['65_group']=='B')&(coord_sorted['227_group']=='B')&(coord_sorted['294_group']=='A')].index
intersect_65_227 = coord_sorted[(coord_sorted['65_group']=='A')&(coord_sorted['227_group']=='A')&(coord_sorted['294_group']=='B')].index
intersect_65_294 = coord_sorted[(coord_sorted['65_group']=='A')&(coord_sorted['227_group']=='B')&(coord_sorted['294_group']=='A')].index
intersect_227_294 = coord_sorted[(coord_sorted['65_group']=='B')&(coord_sorted['227_group']=='A')&(coord_sorted['294_group']=='A')].index
intersect_65_227_294 = coord_sorted[(coord_sorted['65_group']=='A')&(coord_sorted['227_group']=='A')&(coord_sorted['294_group']=='A')].index
#%%
binding_indices=coord_age[(coord_age['65_group']=='A')|(coord_age['227_group']=='A')|(coord_age['294_group']=='A')].index
nonbinding_indices=coord_age[(coord_age['65_group']=='B')&(coord_age['227_group']=='B')&(coord_age['294_group']=='B')].index
phylop_bind = phyloP_matrix[binding_indices]
phylop_nonbind = phyloP_matrix[nonbinding_indices]
alignemnt_bind=filtered_alignment[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=filtered_alignment[nonbinding_indices]
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
plots.plot(
    heatmap=False,
    show_alignment=False,
    alignment=alignment_sorted,
    #logos=True, 
    aggregated=True, 
    aggregated_data=[-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['blue','red'],
    agg_ylim=[[0,10],[-10,0]],
    agg_yhighlight=[[61,69],[224,232],],
    agg_yhighlight_col= ['orange','orange'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_ylabel_right=['conserved\nmotif>no motif','accelerated\nmotif<no motif'], 
    agg_ylabel_right_fs=10,
    agg_ylabel_right_pos=[1.2,0.5],
    agg_ylabel=['-log10(adjusted p-value)',None],
    agg_ylabel_ypos=[0.05,None],
    agg_ylabel_xpos=0,
    agg_ylabel_fs=10,
    #agg_yscale=['log','log'],
    agg_plot_title=['phyloP of THE1C, with both major CLOCK motifs',None],
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
#%%
age_of_interest_list=[0.0,6.7,9.06,15.76,20.19,29.44,43.2,76.0,105.0]
for idx, age_of_interest in enumerate(age_of_interest_list):
    coord_age_filtered=coord_age[coord_age['te_age']==age_of_interest]
    binding_indices = coord_age_filtered[(coord_age_filtered['right_group']=='A')&(coord_age_filtered['left_group']=='A')].index
    nonbinding_indices = coord_age_filtered[(coord_age_filtered['right_group']=='B')&(coord_age_filtered['left_group']=='B')].index
    phylop_bind = phyloP_matrix[binding_indices]
    phylop_nonbind = phyloP_matrix[nonbinding_indices]
    alignemnt_bind=filtered_alignment[binding_indices]
    phylop_bind[alignemnt_bind == 0] = np.nan
    alignemnt_nonbind=filtered_alignment[nonbinding_indices]
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
coord_age_filtered=coord_age[coord_age['te_age']==age_of_interest]
binding_indices = coord_age_filtered[(coord_age_filtered['right_group']=='A')&(coord_age_filtered['left_group']=='A')].index
nonbinding_indices = coord_age_filtered[(coord_age_filtered['right_group']=='B')&(coord_age_filtered['left_group']=='B')].index
phylop_bind = phyloP_matrix[binding_indices]
phylop_nonbind = phyloP_matrix[nonbinding_indices]
alignemnt_bind=filtered_alignment[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=filtered_alignment[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=fdr_control_with_nans(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=fdr_control_with_nans(p_value_less)
#%%
#group_name = age_ref_table[age_ref_table['age']==age_of_interest]['representative'].values[0]
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot(
    heatmap=False,
    show_alignment=False,
    alignment=alignment_sorted,
    #logos=True, 
    aggregated=True, 
    aggregated_data=[-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['blue','red'],
    agg_ylim=[[0,10],[-10,0]],
    agg_yhighlight=[[61,69],[224,232],],
    agg_yhighlight_col= ['orange','orange'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_ylabel_right=['conserved\nmotif>no motif','accelerated\nmotif<no motif'], 
    agg_ylabel_right_fs=10,
    agg_ylabel_right_pos=[1.2,0.5],
    agg_ylabel=['-log10(adjusted p-value)',None],
    agg_ylabel_ypos=[0.05,None],
    agg_ylabel_xpos=0,
    agg_ylabel_fs=10,
    #agg_plot_title=['phyloP of THE1C, grouped by AP-1 motif',None],
    agg_plot_title_fs= 12, 
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    agg_plottext=[None,f'new-world monkey (43.2 MYA)\nTE w/ motifs: {len(binding_indices)}\nTE w/o motifs: {len(nonbinding_indices)}'],
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
