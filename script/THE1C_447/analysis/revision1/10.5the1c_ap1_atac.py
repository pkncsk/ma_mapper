#%%
from curses import meta
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
#%%
atac_df_bed = pd.read_csv('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/phenotype_tracks/atac_coord_hg38.bed', sep='\t', header = None)
atac=mapper.map_and_overlay(alignment_file, coord_file, atac_df_bed, data_format='bed', custom_id=True, strand_overlap=False, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
atac=atac[noNA_indices]
#%%
import importlib
importlib.reload(mapper)
mean_phylop=mapper.normalise(alignment=alignment_filtered, mapped_data=phylop_447)
cov_atac=mapper.normalise(alignment=alignment_filtered, mapped_data=atac, method='perc_coverage')
cov_ap1=mapper.normalise(alignment=alignment_filtered, mapped_data=ap1, method='perc_coverage')
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [atac], 
    alignment=alignment_filtered,
    show_alignment=False, 
    heatmap_color=['Oranges','Greens'],
    heatmap_mode='overlay',
    heatmap_title=['ATAC-seq signal peaks on THE1C MSA'],
    heatmap_title_fs = 10, 
    anno_ylabel = 'sequences',
    anno_ylabel_fs=10,
    vlim = [[0,1]], 
    opacity = 0.9, 
    hm_transparency_mode = 'gradient',
    aggregated=True, 
    aggregated_data=[cov_atac], 
    agg_colset=['orange',],
    agg_ylim=[[0,15]],
    #agg_titles=['ap1 motif'], 
    agg_ylabel=['perc_coverage'],
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
#find peaks
import scipy
peaks, _ = scipy.signal.find_peaks(cov_ap1, width = 6)
#%%
#right,left
#highest_peak_index = peaks[np.argmax(mean_znf808[peaks])]
binding_indices = np.unique(np.where(ap1[:, peaks] != 0)[0])
nonbinding_indices=list(set(np.arange(ap1.shape[0])) - set(binding_indices))
#%%
binding_indices_right = np.unique(np.where(ap1[:, 306] != 0)[0])
nonbinding_indices_right=list(set(np.arange(ap1.shape[0])) - set(binding_indices_right))
binding_indices_left = np.unique(np.where(ap1[:, 100] != 0)[0])
nonbinding_indices_left=list(set(np.arange(ap1.shape[0])) - set(binding_indices_left))
#%%
metadata_age['right_group'] = 'B'
metadata_age.loc[metadata_age.index.isin(binding_indices_right), 'right_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices_right), 'right_group'] = 'B'
metadata_age['left_group'] = 'B'
metadata_age.loc[metadata_age.index.isin(binding_indices_left), 'left_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices_left), 'left_group'] = 'B'
#%%
metadata_age['motif_group'] = 'D'
metadata_age.loc[(metadata_age['right_group'] == 'A') & (metadata_age['left_group'] == 'A'), 'motif_group'] = 'A'
metadata_age.loc[(metadata_age['right_group'] == 'A') & (metadata_age['left_group'] == 'B'), 'motif_group'] = 'B'
metadata_age.loc[(metadata_age['right_group'] == 'B') & (metadata_age['left_group'] == 'A'), 'motif_group'] = 'C'
metadata_age.loc[(metadata_age['right_group'] == 'B') & (metadata_age['left_group'] == 'B'), 'motif_group'] = 'D'

#%%
import scipy
peaks, _ = scipy.signal.find_peaks(cov_atac, width = 6)
# %%
annot_indices = np.unique(np.where(atac[:, peaks] != 0)[0])
nonannot_indices=list(set(np.arange(atac.shape[0])) - set(annot_indices))
#%%
metadata_age['atac_group'] = 'B'
metadata_age.loc[metadata_age.index.isin(annot_indices), 'atac_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(nonannot_indices), 'atac_group'] = 'B'
atac_sort=metadata_age.sort_values(['atac_group']).index


#%%
test_set=metadata_age.sort_values(['atac_group']).reset_index()
test_set[(test_set['atac_group']==1)].index
# %%
import pybedtools
atac_df_bed = pd.read_csv('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/phenotype_tracks/atac_coord_hg38.bed', sep='\t', header = None)
atac_bed=pybedtools.BedTool.from_dataframe(atac_df_bed)
coord_bed=pybedtools.BedTool.from_dataframe(coord_file)
intersect_bed=coord_bed.intersect(atac_bed,loj=True, wa = True, wb =True, s=False)
intersect_df=intersect_bed.to_dataframe(names=['chrom', 'start', 'end', 'name', 'score', 'strand','meta_id','chrom2', 'start2', 'end2', 'name_peak', 'score2', 'strand2',])
intersect_df=intersect_df[['chrom', 'start', 'end', 'name', 'score', 'strand','name_peak']]
# %%
atac_table = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/phenotype_tracks/encode_v3_score_matrix_norm_class_kzfps_embed_all_clust_TSS.csv'
atac_df=pd.read_csv(atac_table, sep=',', low_memory=False, index_col = 0)
atac_the1_table = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/phenotype_tracks/encode_v3_peaks_the1_counts_sort_kzfps_sort.csv'
atac_the1_df=pd.read_csv(atac_the1_table, sep=',', low_memory=False)
#%%
atac_the1_df_lean = atac_the1_df[['name','cluster_id','TSS','TSS_5000']].rename(columns={'name':'name_peak'})
# %%
intersect_info=intersect_df.merge(atac_the1_df_lean, on='name_peak')
# %%
metadata_info=metadata_age.merge(intersect_info, on=['chrom','start','end','name','score','strand'], how='left')
#%%
metadata_info['immune_cluster'] = 'B'
metadata_info.loc[metadata_info['cluster_id'].isin([2,3,12,14,23,34,]), 'immune_cluster'] = 'A'

# %%
metadata_sorted=metadata_info.sort_values(['atac_group','te_age'])
metadata_sorted['immune_cluster'] = 'B'
metadata_sorted.loc[metadata_sorted['cluster_id'].isin([2,3,12,14,23,34,]), 'immune_cluster'] = 'A'
sorted_indices=metadata_sorted.index
age_subgroups = np.unique(metadata_sorted['te_age'].sort_values())
age_subgroup = {subgroup: num for num, subgroup in enumerate(age_subgroups)}
age_anno=metadata_sorted['te_age'].map(age_subgroup)
motif_anno = metadata_sorted['motif_group']
atac_anno =metadata_sorted['atac_group']
tss_anno = metadata_sorted['TSS']
tss5k_anno = metadata_sorted['TSS_5000']
imm_anno = metadata_sorted['immune_cluster']
phylop_sorted=phylop_447[sorted_indices]
alignment_sorted=alignment_filtered[sorted_indices]
atac_sorted = atac[sorted_indices]
#%%
test_set=metadata_sorted.reset_index()
test_set[(test_set['atac_group']=='A')].index
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [atac_sorted], 
    alignment=alignment_sorted,
    show_alignment=True, 
    heatmap_color=['Oranges'],
    heatmap_mode='overlay', 
    vlim = [[0,1]], 
    opacity = 0.8, 
        hm_transparency_mode = 'gradient',
    heatmap_title=['ATAC-seq signal peaks on THE1C, grouped by ATAC-seq signals'],
    heatmap_title_fs = 10, 
    heatmap_xhighlight=[[1278,1278]],
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
    hm_interpolation=None,
    hm_transparency_mode='gradient', 
    annotation=True, 
    anno_col = ['Blues'], 
    annotation_data=[age_anno],
    anno_cbar_label=[age_subgroups],
    anno_cbar_title=['TEA-TIME'], 
    anno_cbar_even_pos=-0.19,
    colorbar=True,
    colorbar_steps=[0.1],
    heatmap_title=['phyloP of THE1C, grouped by ATAC-seq signals'],
    heatmap_title_fs = 10, 
    heatmap_xlabel = 'position (bp)',
    heatmap_xlabel_fs=10,
    anno_ylabel = 'sequences',
    anno_ylabel_fs=10,
    heatmap_xhighlight=[[1278,1278]],
    heatmap_xhighlight_col= ['black','black','black'],
    heatmap_xhighlight_alpha=[0.8,0.8,0.8],
    )
#%%
binding_indices=metadata_sorted[(metadata_sorted['atac_group']=='A')].index
nonbinding_indices=metadata_sorted[(metadata_sorted['atac_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
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
    #agg_yhighlight=[[100,350],[212,232]],
    #agg_yhighlight_col= ['purple','purple'],
    #agg_yhighlight_alpha=[0.1,0.1],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['conserved\nsignal>no signal','accelerated\nsignal<no signal'], 
    agg_titles_fs=10,
    agg_titles_pos=[1.2,0.5],
    agg_ylabel=['-log10(adjusted p-value)',None],
    agg_ylabel_ypos=[0.05,None],
    agg_ylabel_xpos=0,
    agg_ylabel_fs=10,
    agg_plot_title=['phyloP of THE1C, grouped by ATAC-seq signals',None],
    agg_plot_title_fs= 12, 
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    agg_plottext=[None, f'TE w/ signals: {len(binding_indices)}\nTE w/o signals: {len(nonbinding_indices)}'],
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
   agg_yhighlight=[[100,350],[212,232],[98,108],[304,314],[109,122],[175,190],[63,68],[224,232],],
    agg_yhighlight_col= ['purple','purple','green','green','blue','blue','red','red'],
    agg_yhighlight_alpha=[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['conserved\nsignal>no signal','accelerated\nsignal<no signal'], 
    agg_titles_fs=10,
    agg_titles_pos=[1.2,0.5],
    agg_ylabel=['-log10(adjusted p-value)',None],
    agg_ylabel_ypos=[0.05,None],
    agg_ylabel_xpos=0,
    agg_ylabel_fs=10,
    agg_plot_title=['phyloP of THE1C, grouped by ATAC-seq signals',None],
    agg_plot_title_fs= 12, 
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    agg_plottext=[None, f'TE w/ signals: {len(binding_indices)}\nTE w/o signals: {len(nonbinding_indices)}'],
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
age_of_interest_list=config.age_canon[0:7]
for idx, age_of_interest in enumerate(age_of_interest_list):
    metadata_age_filtered=metadata_sorted[metadata_sorted['te_age']==age_of_interest]
    binding_indices=metadata_age_filtered[(metadata_age_filtered['atac_group']==1)].index
    nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['atac_group']==0)].index
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
age_of_interest=43.2
metadata_age_filtered=metadata_sorted[metadata_sorted['te_age']==age_of_interest]
binding_indices=metadata_age_filtered[(metadata_age_filtered['atac_group']==1)].index
nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['atac_group']==0)].index
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
   # agg_yhighlight=[[100,350],[212,232]],
    #agg_yhighlight_col= ['purple','purple'],
    #agg_yhighlight_alpha=[0.1,0.1,],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
   agg_titles=['conserved\nsignal>no signal','accelerated\nsignal<no signal'], 
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
    agg_plottext=[None, f'{group_name} ({age_of_interest} MYA)\nTE w/ signals: {len(binding_indices)}\nTE w/o signals: {len(nonbinding_indices)}'],
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
metadata_sorted=metadata_sorted.sort_values(['right_group','immune_cluster','te_age'])
sorted_indices=metadata_sorted.index
age_subgroups = np.unique(metadata_sorted['te_age'].sort_values())
age_subgroup = {subgroup: num for num, subgroup in enumerate(age_subgroups)}
age_anno=metadata_sorted['te_age'].map(age_subgroup)
motif_anno = metadata_sorted['right_group']
atac_anno =metadata_sorted['atac_group']
tss_anno = metadata_sorted['TSS']
tss5k_anno = metadata_sorted['TSS_5000']
imm_anno = metadata_sorted['immune_cluster']
phylop_sorted=phylop_447[sorted_indices]
alignment_sorted=alignment_filtered[sorted_indices]
atac_sorted = atac[sorted_indices]
ap1_sorted = ap1[sorted_indices]
#%%
metadata_sorted['motif_group2'] = 'D'
metadata_sorted.loc[(metadata_sorted['right_group'] == 'A') & (metadata_sorted['immune_cluster'] == 'A'), 'motif_group2'] = 'A'
metadata_sorted.loc[(metadata_sorted['right_group'] == 'A') & (metadata_sorted['immune_cluster'] == 'B'), 'motif_group2'] = 'B'
metadata_sorted.loc[(metadata_sorted['right_group'] == 'B') & (metadata_sorted['immune_cluster'] == 'A'), 'motif_group2'] = 'C'
metadata_sorted.loc[(metadata_sorted['right_group'] == 'B') & (metadata_sorted['immune_cluster'] == 'B'), 'motif_group2'] = 'D'
subgroups = np.unique(metadata_sorted['motif_group2'].astype(str))
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
motif_anno2=metadata_sorted['motif_group2'].map(numerical_subgroup)
#%%
test_set=metadata_sorted.reset_index()
test_set[(test_set['right_group']=='A')&(test_set['immune_cluster']=='A')].index
#%%
test_set[(test_set['right_group']=='A')&(test_set['immune_cluster']=='B')].index
#%%
test_set[(test_set['right_group']=='B')&(test_set['immune_cluster']=='A')].index
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [atac_sorted, ap1_sorted], 
    alignment=alignment_sorted,
    show_alignment=True, 
    heatmap_color=['Oranges', 'Greens'],
    heatmap_mode='overlay', 
    vlim = [[0,1]], 
    opacity = 0.8, 
        hm_transparency_mode = 'gradient',
    heatmap_title=['ATAC-seq signals and AP-1 motifs on THE1C,\ngrouped by both motifs and signals'],
    heatmap_title_fs = 10, 
    heatmap_xhighlight=[[361,361],[5573,5573],[5713,5713]],
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
    hm_interpolation = None, 
     hm_transparency_mode = 'gradient',
    annotation=True, 
    anno_col = ['Blues',['purple','red','blue','white']], 
    annotation_data=[age_anno,motif_anno2],
    anno_cbar_label=[age_subgroups,['AP-1 and\naccessible','only AP-1','only accessible','no signals\ninaccessible']],
    anno_cbar_title=['TEA-TIME','AP-1 and\naccessibility\n'], 
    anno_cbar_even_pos=-0.3,
    colorbar=True,
    colorbar_steps=[0.1],
    heatmap_title=['phyloP of THE1C, grouped by both motifs and signals'],
    heatmap_title_fs = 10, 
    heatmap_xlabel = 'position (bp)',
    heatmap_xlabel_fs=10,
    anno_ylabel = 'sequences',
    anno_ylabel_fs=10,
      heatmap_xhighlight=[[361,361],[5573,5573],[5713,5713]],
    heatmap_xhighlight_col= ['black','black','black'],
    heatmap_xhighlight_alpha=[0.8,0.8,0.8],
    )
# %%
#%%
binding_indices=metadata_sorted[(metadata_sorted['right_group']=='A')&(metadata_sorted['immune_cluster']=='A')].index
nonbinding_indices=metadata_sorted[(metadata_sorted['right_group']=='B')&(metadata_sorted['immune_cluster']=='B')].index
phylop_bind = phylop_447[binding_indices]
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
    agg_yhighlight=[[304,314]],
    agg_yhighlight_col= ['green','purple'],
    agg_yhighlight_alpha=[0.1,0.1],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['conserved\nsignal>no signal','accelerated\nsignal<no signal'], 
    agg_titles_fs=10,
    agg_titles_pos=[1.2,0.5],
    agg_ylabel=['-log10(adjusted p-value)',None],
    agg_ylabel_ypos=[0.05,None],
    agg_ylabel_xpos=0,
    agg_ylabel_fs=10,
    agg_plot_title=['phyloP of THE1C, with AP-1 motifs and accessibility in immune cells',None],
    agg_plot_title_fs= 12, 
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    agg_plottext=[None, f'TE w/ signals: {len(binding_indices)}\nTE w/o signals: {len(nonbinding_indices)}'],
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
from ma_mapper import sequence_alignment
maybe_motif=alignemnt_bind[:,300:320]
seq=sequence_alignment.parsed_array_to_sequence(maybe_motif, prefix='THE1C', output_filepath='/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/THE1C_447/analysis/revision1/motif.fasta')
# %%
age_of_interest_list=config.age_canon[0:7]
for idx, age_of_interest in enumerate(age_of_interest_list):
    metadata_age_filtered=metadata_sorted[metadata_sorted['te_age']==age_of_interest]
    binding_indices=metadata_age_filtered[(metadata_age_filtered['right_group']==1)&(metadata_age_filtered['immune_cluster']==1)].index
    nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['right_group']==0)&(metadata_age_filtered['immune_cluster']==0)].index
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
age_of_interest=43.2
metadata_age_filtered=metadata_sorted[metadata_sorted['te_age']==age_of_interest]
binding_indices=metadata_age_filtered[(metadata_age_filtered['right_group']==1)&(metadata_age_filtered['immune_cluster']==1)].index
nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['right_group']==0)&(metadata_age_filtered['immune_cluster']==0)].index
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
    agg_yhighlight=[[304,314]],
    agg_yhighlight_col= ['green','purple'],
    agg_yhighlight_alpha=[0.1,0.1,],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
   agg_titles=['conserved\nsignal>no signal','accelerated\nsignal<no signal'], 
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
    agg_plottext=[None, f'{group_name} ({age_of_interest} MYA)\nTE w/ signals: {len(binding_indices)}\nTE w/o signals: {len(nonbinding_indices)}'],
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
#%%
binding_indices=metadata_sorted[(metadata_sorted['right_group']=='A')&(metadata_sorted['immune_cluster']=='B')].index
nonbinding_indices=metadata_sorted[(metadata_sorted['right_group']=='B')&(metadata_sorted['immune_cluster']=='B')].index
phylop_bind = phylop_447[binding_indices]
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
    agg_ylim=[[0,25],[-25,0]],
    agg_yhighlight=[[304,314]],
    agg_yhighlight_col= ['green','purple'],
    agg_yhighlight_alpha=[0.1,0.1],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['conserved\nsignal>no signal','accelerated\nsignal<no signal'], 
    agg_titles_fs=10,
    agg_titles_pos=[1.2,0.5],
    agg_ylabel=['-log10(adjusted p-value)',None],
    agg_ylabel_ypos=[0.05,None],
    agg_ylabel_xpos=0,
    agg_ylabel_fs=10,
    agg_plot_title=['phyloP of THE1C, with only AP-1 motifs',None],
    agg_plot_title_fs= 12, 
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    agg_plottext=[None, f'TE w/ signals: {len(binding_indices)}\nTE w/o signals: {len(nonbinding_indices)}'],
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
# %%
age_of_interest_list=config.age_canon[0:7]
for idx, age_of_interest in enumerate(age_of_interest_list):
    metadata_age_filtered=metadata_sorted[metadata_sorted['te_age']==age_of_interest]
    binding_indices=metadata_age_filtered[(metadata_age_filtered['right_group']==1)&(metadata_age_filtered['immune_cluster']==0)].index
    nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['right_group']==0)&(metadata_age_filtered['immune_cluster']==0)].index
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
age_of_interest=43.2
metadata_age_filtered=metadata_sorted[metadata_sorted['te_age']==age_of_interest]
binding_indices=metadata_age_filtered[(metadata_age_filtered['right_group']==1)&(metadata_age_filtered['immune_cluster']==0)].index
nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['right_group']==0)&(metadata_age_filtered['immune_cluster']==0)].index
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
    agg_ylim=[[0,25],[-25,0]],
    agg_yhighlight=[[304,314]],
    agg_yhighlight_col= ['green','purple'],
    agg_yhighlight_alpha=[0.1,0.1,],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
   agg_titles=['conserved\nsignal>no signal','accelerated\nsignal<no signal'], 
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
    agg_plottext=[None, f'{group_name} ({age_of_interest} MYA)\nTE w/ signals: {len(binding_indices)}\nTE w/o signals: {len(nonbinding_indices)}'],
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
binding_indices=metadata_sorted[(metadata_sorted['right_group']==0)&(metadata_sorted['immune_cluster']==1)].index
nonbinding_indices=metadata_sorted[(metadata_sorted['right_group']==0)&(metadata_sorted['immune_cluster']==0)].index
phylop_bind = phylop_447[binding_indices]
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
    agg_ylim=[[0,5],[-5,0]],
    agg_yhighlight=[[304,314]],
    agg_yhighlight_col= ['green','purple'],
    agg_yhighlight_alpha=[0.1,0.1],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['conserved\nsignal>no signal','accelerated\nsignal<no signal'], 
    agg_titles_fs=10,
    agg_titles_pos=[1.2,0.5],
    agg_ylabel=['-log10(adjusted p-value)',None],
    agg_ylabel_ypos=[0.05,None],
    agg_ylabel_xpos=0,
    agg_ylabel_fs=10,
    agg_plot_title=['phyloP of THE1C, with AP-1 motifs and accessibility in immune cells',None],
    agg_plot_title_fs= 12, 
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    agg_plottext=[None, f'TE w/ signals: {len(binding_indices)}\nTE w/o signals: {len(nonbinding_indices)}'],
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
# %%
age_of_interest_list=config.age_canon[0:7]
for idx, age_of_interest in enumerate(age_of_interest_list):
    metadata_age_filtered=metadata_sorted[metadata_sorted['te_age']==age_of_interest]
    binding_indices=metadata_age_filtered[(metadata_age_filtered['right_group']==0)&(metadata_age_filtered['immune_cluster']==1)].index
    nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['right_group']==0)&(metadata_age_filtered['immune_cluster']==0)].index
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
age_of_interest=43.2
metadata_age_filtered=metadata_sorted[metadata_sorted['te_age']==age_of_interest]
binding_indices=metadata_age_filtered[(metadata_age_filtered['right_group']==0)&(metadata_age_filtered['immune_cluster']==1)].index
nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['right_group']==0)&(metadata_age_filtered['immune_cluster']==0)].index
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
    agg_ylim=[[0,25],[-25,0]],
    agg_yhighlight=[[304,314]],
    agg_yhighlight_col= ['green','purple'],
    agg_yhighlight_alpha=[0.1,0.1,],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
   agg_titles=['conserved\nsignal>no signal','accelerated\nsignal<no signal'], 
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
    agg_plottext=[None, f'{group_name} ({age_of_interest} MYA)\nTE w/ signals: {len(binding_indices)}\nTE w/o signals: {len(nonbinding_indices)}'],
    agg_plottext_fs=10,
    agg_plottext_pos=[0.99,0.01],
    agg_major_tick=20,
    colorbar=False,
    agg_h=20,
    figsize= [80,40],
    #xlim=[304,313],
    #gg_major_tick=1,
    )
