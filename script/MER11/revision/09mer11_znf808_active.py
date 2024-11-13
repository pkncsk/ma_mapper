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
subfamily='MER11'
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
#%%
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file,col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
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
metadata_age['subfam'] = metadata_age['group'].str.split('_').str[0]
subfam_subgroups = np.unique(metadata_age['subfam'].astype(str))
subfam_subgroup = {subgroup: num for num, subgroup in enumerate(subfam_subgroups)}
subfam_anno=metadata_age['subfam'].map(subfam_subgroup)
#%%
coord_file=mapper.extract_metadata_from_alignment(alignment_file)
#%%
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/zoonomia447/hg38.phyloP447way.bw'
phylop_447=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
phylop_447=phylop_447[noNA_indices]
# %%
znf808active_table = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/phenotype_tracks/znf808_active_hg38.bed'
znf808a_df=pd.read_csv(znf808active_table, sep='\t', low_memory=False)
znf808a_df_bed = znf808a_df.iloc[:,0:3]
znf808a_df_bed.columns = ['chrom','start','end']
znf808a_df_bed['name'] = 'ZNF808a_' + znf808a_df_bed.index.astype(str)
znf808a_df_bed['score']=10
znf808a_df_bed['strand'] = '.'
znf808a=mapper.map_and_overlay(alignment_file, coord_file, znf808a_df_bed, data_format='bed', custom_id=True, strand_overlap = False, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
znf808a=znf808a[noNA_indices]
#%%
import importlib
importlib.reload(mapper)
mean_phylop=mapper.normalise(alignment=alignment_filtered, mapped_data=phylop_447)
cov_znf808a=mapper.normalise(alignment=alignment_filtered, mapped_data=znf808a, method='perc_coverage')
#%%
# %%
#find peaks
import scipy
peaks, _ = scipy.signal.find_peaks(cov_znf808a, width = 6)
#%%
#highest_peak_index = peaks[np.argmax(mean_znf808[peaks])]
binding_indices = np.unique(np.where(znf808a[:, peaks] != 0)[0])
nonbinding_indices=list(set(np.arange(znf808a.shape[0])) - set(binding_indices))
metadata_age['znf808a_group'] = 'B'
metadata_age.loc[metadata_age.index.isin(binding_indices), 'znf808a_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices), 'znf808a_group'] = 'B'
#%%
import pybedtools
znf808a_bed=pybedtools.BedTool.from_dataframe(znf808a_df_bed)
coord_bed=pybedtools.BedTool.from_dataframe(coord_file)
intersect_bed=coord_bed.intersect(znf808a_bed,loj=True, c=True, s=False)
# %%
intersect_df=intersect_bed.to_dataframe(names=['chrom', 'start', 'end', 'name', 'score', 'strand','meta_id','znf808a',])
intersect_df=intersect_df[['chrom', 'start', 'end', 'name', 'score', 'strand','znf808a']]
# %%
metadata_age_lean=metadata_age.merge(intersect_df, on=['chrom', 'start', 'end', 'name', 'score', 'strand'], how = 'left')
# %%
#metadata_sorted=metadata_age_lean.sort_values(['subfam'])
metadata_sorted=metadata_age_lean.sort_values(['subfam','znf808a_group','te_age'])
sorted_indices = metadata_sorted.index
znf808a_sorted = znf808a[sorted_indices]
subfam_anno_sorted = subfam_anno[sorted_indices]
age_anno_sorted = age_anno[sorted_indices]
znf808a_anno_sorted = metadata_sorted['znf808a_group']
znf808a_anno = metadata_sorted['znf808a']
phylop_sorted = phylop_447[sorted_indices]
alignment_sorted = alignment_filtered[sorted_indices]
#%%
test_set=metadata_sorted.reset_index()
#%%
test_set[(test_set['subfam']=='MER11A')].index
#%%
test_set[(test_set['subfam']=='MER11B')].index
#%%
test_set[(test_set['subfam']=='MER11C')].index
#%%
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [znf808a_sorted], 
    alignment=alignment_sorted,
    show_alignment=True, 
    heatmap_color=['Greens','Blues'],
    heatmap_mode='overlay',
    heatmap_title=['H3K27ac signal peaks on MER11 MSA, sorted by H3K27ac signals'],
    heatmap_title_fs = 10, 
    anno_ylabel = 'sequences',
    anno_ylabel_fs=10,
    vlim = [[0,7]], 
    opacity = 0.5, 
    hm_transparency_mode = 'gradient',
    aggregated=False, 
    aggregated_data=[cov_znf808a], 
    agg_colset=['green'],
    agg_ylim=[[0,15],[0,110]],
    agg_titles=['total\ncoverage','coverage\nper subfamily'], 
        agg_titles_fs=8,
    agg_titles_pos=[1.2,0.5],
    agg_ylabel=['perc_coverage',None],
    agg_ylabel_ypos=[0.05,None],
     agg_ylabel_xpos=0,
    agg_ylabel_fs=8,
    agg_xlabel='position (bp)',
     agg_xlabel_fs=10,
    agg_major_tick=100,
        annotation=False, 
    anno_col = ['Blues',['red', 'yellow', 'blue']], 
    #anno_title=['TEA-TIME','subfamily'],
    annotation_data=[age_anno_sorted,subfam_anno_sorted],
    anno_cbar_label=[age_subgroups,['MER11A','MER11B','MER11C']],
    anno_cbar_title=['TEA-TIME','subfamily'],  
    anno_cbar_even_pos=-0.19,
        heatmap_xhighlight=[[940,940],[1475,1475]],
    heatmap_xhighlight_col= ['black','black','black'],
    heatmap_xhighlight_alpha=[0.5,0.5,0.5],
    heatmap_xlabel='position (bp)',
heatmap_xlabel_fs=10,
    )
#%%
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
    anno_col = ['Blues',['red', 'yellow', 'blue']], 
    #anno_title=['TEA-TIME','signals','subfamily'],
    annotation_data=[age_anno_sorted,subfam_anno_sorted],
    anno_cbar_label=[age_subgroups,['MER11A','MER11B','MER11C']],
    anno_cbar_title=['TEA-TIME','signals','subfamily'],  
    anno_cbar_even_pos=-0.19,
    colorbar=True,
    colorbar_steps=[0.1],
    heatmap_title=['phyloP of THE1C, grouped by H3K27ac signals'],
    heatmap_title_fs = 10, 
    heatmap_xlabel = 'position (bp)',
    heatmap_xlabel_fs=10,
    anno_ylabel = 'sequences',
    anno_ylabel_fs=10,
     heatmap_xhighlight=[[940,940],[1475,1475]],
    heatmap_xhighlight_col= ['black','black','black'],
    heatmap_xhighlight_alpha=[0.8,0.8,0.8],
    agg_major_tick=100,
    )
#%%
subfamily_of_interest = 'MER11A'
metadata_age_filtered=metadata_age[metadata_age['subfam']==subfamily_of_interest]
binding_indices=metadata_age_filtered[(metadata_age_filtered['znf808a_group']=='A')].index
nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['znf808a_group']=='B')].index
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
    agg_ylim=[[0,5],[-5,0]],
    agg_yhighlight=[[450,800],[600,625]],
    agg_yhighlight_col= ['green','green'],
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
    agg_plot_title=[f'phyloP of {subfamily_of_interest}, grouped by H3K27ac signals',None],
    agg_plot_title_fs= 12, 
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    agg_plottext=[None, f'TE w/ signals: {len(binding_indices)}\nTE w/o signals: {len(nonbinding_indices)}'],
    agg_plottext_fs=10,
    agg_plottext_pos=[0.99,0.01],
    agg_major_tick=100,
    colorbar=False,
    agg_h=20,
    figsize= [80,40],
    #xlim=[304,313],
    #gg_major_tick=1,
    )
#%%
age_of_interest_list=config.age_canon[0:6]
subfamily_of_interest = 'MER11A'
for idx, age_of_interest in enumerate(age_of_interest_list):
    metadata_age_filtered_=metadata_age[metadata_age['subfam']==subfamily_of_interest]
    metadata_age_filtered=metadata_age_filtered_[metadata_age_filtered_['te_age']==age_of_interest]
    binding_indices = metadata_age_filtered[(metadata_age_filtered['znf808a_group']=='A')].index
    nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['znf808a_group']=='B')].index
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
#%%
subfamily_of_interest = 'MER11B'
metadata_age_filtered=metadata_age[metadata_age['subfam']==subfamily_of_interest]
binding_indices=metadata_age_filtered[(metadata_age_filtered['znf808a_group']=='A')].index
nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['znf808a_group']=='B')].index
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
    agg_yhighlight=[[450,800],[600,625]],
    agg_yhighlight_col= ['green','green'],
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
    agg_plot_title=[f'phyloP of {subfamily_of_interest}, grouped by H3K27ac signals',None],
    agg_plot_title_fs= 12, 
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    agg_plottext=[None, f'TE w/ signals: {len(binding_indices)}\nTE w/o signals: {len(nonbinding_indices)}'],
    agg_plottext_fs=10,
    agg_plottext_pos=[0.99,0.01],
    agg_major_tick=100,
    colorbar=False,
    agg_h=20,
    figsize= [80,40],
    #xlim=[304,313],
    #gg_major_tick=1,
    )
#%%
age_of_interest_list=config.age_canon[0:6]
subfamily_of_interest = 'MER11B'
for idx, age_of_interest in enumerate(age_of_interest_list):
    metadata_age_filtered_=metadata_age[metadata_age['subfam']==subfamily_of_interest]
    metadata_age_filtered=metadata_age_filtered_[metadata_age_filtered_['te_age']==age_of_interest]
    binding_indices = metadata_age_filtered[(metadata_age_filtered['znf808a_group']=='A')].index
    nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['znf808a_group']=='B')].index
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
# %%
age_of_interest=20.19
subfamily_of_interest = 'MER11B'
metadata_age_filtered_=metadata_age[metadata_age['subfam']==subfamily_of_interest]
metadata_age_filtered=metadata_age_filtered_[metadata_age_filtered_['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['znf808a_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['znf808a_group']=='B')].index
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
    agg_yhighlight=[[450,800],[600,625]],
    agg_yhighlight_col= ['green','green'],
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
    agg_major_tick=100,
    colorbar=False,
    agg_h=20,
    figsize= [80,40],
    #xlim=[304,313],
    #gg_major_tick=1,
    )
# %%
#%%
subfamily_of_interest = 'MER11C'
metadata_age_filtered=metadata_age[metadata_age['subfam']==subfamily_of_interest]
binding_indices=metadata_age_filtered[(metadata_age_filtered['znf808a_group']=='A')].index
nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['znf808a_group']=='B')].index
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
    agg_yhighlight=[[450,800],[600,625]],
    agg_yhighlight_col= ['green','green'],
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
    agg_plot_title=[f'phyloP of {subfamily_of_interest}, grouped by H3K27ac signals',None],
    agg_plot_title_fs= 12, 
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    agg_plottext=[None, f'TE w/ signals: {len(binding_indices)}\nTE w/o signals: {len(nonbinding_indices)}'],
    agg_plottext_fs=10,
    agg_plottext_pos=[0.99,0.01],
    agg_major_tick=100,
    colorbar=False,
    agg_h=20,
    figsize= [80,40],
    #xlim=[304,313],
    #gg_major_tick=1,
    )
#%%
age_of_interest_list=config.age_canon[0:6]
subfamily_of_interest = 'MER11C'
for idx, age_of_interest in enumerate(age_of_interest_list):
    metadata_age_filtered_=metadata_age[metadata_age['subfam']==subfamily_of_interest]
    metadata_age_filtered=metadata_age_filtered_[metadata_age_filtered_['te_age']==age_of_interest]
    binding_indices = metadata_age_filtered[(metadata_age_filtered['znf808a_group']=='A')].index
    nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['znf808a_group']=='B')].index
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
# %%
age_of_interest=20.19
subfamily_of_interest = 'MER11C'
metadata_age_filtered_=metadata_age[metadata_age['subfam']==subfamily_of_interest]
metadata_age_filtered=metadata_age_filtered_[metadata_age_filtered_['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['znf808a_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['znf808a_group']=='B')].index
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
    agg_yhighlight=[[450,800],[600,625]],
    agg_yhighlight_col= ['green','green'],
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
    agg_major_tick=100,
    colorbar=False,
    agg_h=20,
    figsize= [80,40],
    #xlim=[304,313],
    #gg_major_tick=1,
    )
# %%
