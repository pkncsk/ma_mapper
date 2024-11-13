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
bed_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/kzfp_peak_bed/hg38_kzfps_combined.bed'
kzfp_df=pd.read_csv(bed_filepath, sep='\t', header=None)
kzfp_df.columns=['chrom','start','end','name','score','strand']
bed_file=kzfp_df[kzfp_df['name'].str.contains('ZNF808')]
znf808=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap = False, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
znf808=znf808[noNA_indices]
bed_file=kzfp_df[kzfp_df['name'].str.contains('ZNF525')]
znf525=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap = False, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
znf525=znf525[noNA_indices]
#%%
import importlib
importlib.reload(mapper)
mean_phylop=mapper.normalise(alignment=alignment_filtered, mapped_data=phylop_447)
cov_znf808=mapper.normalise(alignment=alignment_filtered, mapped_data=znf808, method='perc_coverage')
cov_znf525=mapper.normalise(alignment=alignment_filtered, mapped_data=znf525, method='perc_coverage')
#%%
# %%
#find peaks
import scipy
peaks, _ = scipy.signal.find_peaks(cov_znf808, width = 6)
#%%
#highest_peak_index = peaks[np.argmax(mean_znf808[peaks])]
binding_indices = np.unique(np.where(znf808[:, peaks] != 0)[0])
nonbinding_indices=list(set(np.arange(znf808.shape[0])) - set(binding_indices))
metadata_age['znf808_group'] = 'No Group'
metadata_age.loc[metadata_age.index.isin(binding_indices), 'znf808_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices), 'znf808_group'] = 'B'
znf808_subgroups = np.unique(metadata_age['znf808_group'].astype(str))
znf808_subgroup = {subgroup: num for num, subgroup in enumerate(znf808_subgroups)}
znf808_subgroup_anno=metadata_age['znf808_group'].map(znf808_subgroup)

# %%
#find peaks
import scipy
peaks, _ = scipy.signal.find_peaks(cov_znf525, width = 6)
#%%
#highest_peak_index = peaks[np.argmax(mean_znf525[peaks])]
binding_indices = np.unique(np.where(znf525[:, peaks] != 0)[0])
nonbinding_indices=list(set(np.arange(znf525.shape[0])) - set(binding_indices))
metadata_age['znf525_group'] = 'No Group'
metadata_age.loc[metadata_age.index.isin(binding_indices), 'znf525_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices), 'znf525_group'] = 'B'
znf525_subgroups = np.unique(metadata_age['znf525_group'].astype(str))
znf525_subgroup = {subgroup: num for num, subgroup in enumerate(znf525_subgroups)}
znf525_subgroup_anno=metadata_age['znf525_group'].map(znf525_subgroup)
#%%
binding_indices=metadata_age[(metadata_age['znf808_group']=='A')|(metadata_age['znf525_group']=='A')].index
nonbinding_indices=metadata_age[(metadata_age['znf808_group']=='B')&(metadata_age['znf525_group']=='B')].index
only_right = metadata_age[(metadata_age['znf808_group']=='A')&(metadata_age['znf525_group']=='B')].index
only_left = metadata_age[(metadata_age['znf808_group']=='B')&(metadata_age['znf525_group']=='A')].index

intersect_right_left = metadata_age[(metadata_age['znf808_group']=='A')&(metadata_age['znf525_group']=='A')].index

#%%
metadata_age['motif_group'] = 'No Group'
metadata_age.loc[metadata_age.index.isin(intersect_right_left), 'motif_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(only_right), 'motif_group'] = 'B'
metadata_age.loc[metadata_age.index.isin(only_left), 'motif_group'] = 'C'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices), 'motif_group'] = 'D'
subgroups = np.unique(metadata_age['motif_group'].astype(str))
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno=metadata_age['motif_group'].map(numerical_subgroup)
#%%
metadata_age_sorted=metadata_age.sort_values(['subfam'])
sorted_indices= metadata_age_sorted.index
znf808_sorted=znf808[sorted_indices]
znf525_sorted=znf525[sorted_indices]
alignment_sorted = alignment_filtered[sorted_indices]
age_anno_sorted=age_anno[sorted_indices]
subfam_anno_sorted = subfam_anno[sorted_indices]
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    alignment=alignment_sorted,
    show_alignment=True, 
    alignment_col='dna', 
    show_alignment_colbar=True,
    heatmap_mode='overlay', 
    opacity = 0.5,
    annotation=True, 
    anno_col = ['Blues',['red', 'yellow', 'blue']], 
    #anno_title=['TEA-TIME','subfamily'],
    annotation_data=[age_anno_sorted,subfam_anno_sorted],
    anno_cbar_label=[age_subgroups,['MER11A','MER11B','MER11C']],
    anno_cbar_title=['TEA-TIME','subfamily'],  
    anno_cbar_even_pos=-0.19,
    colorbar=True,
    colorbar_steps=[0.1],
    heatmap_title=['MER11 alignment, grouped by subfamily'],
    heatmap_title_fs = 10, 
    heatmap_xlabel = 'position (bp)',
    heatmap_xlabel_fs=10,
    anno_ylabel = 'sequences',
    anno_ylabel_fs=10,
    #  heatmap_xhighlight=[[666,666]],
    #heatmap_xhighlight_col= ['black','black','black'],
    #heatmap_xhighlight_alpha=[0.8,0.8,0.8],
    agg_major_tick=100,
    )
#%%
subfamily_of_interest = 'MER11A'
binding_indices=metadata_age[metadata_age['subfam']==subfamily_of_interest].index
znf808_sub = znf808[binding_indices]
alignemnt_sub = alignment_filtered[binding_indices]
cov_znf808_mer11a=mapper.normalise(alignment=alignemnt_sub, mapped_data=znf808_sub, method='perc_coverage')
subfamily_of_interest = 'MER11B'
binding_indices=metadata_age[metadata_age['subfam']==subfamily_of_interest].index
znf808_sub = znf808[binding_indices]
alignemnt_sub = alignment_filtered[binding_indices]
cov_znf808_mer11b=mapper.normalise(alignment=alignemnt_sub, mapped_data=znf808_sub, method='perc_coverage')
subfamily_of_interest = 'MER11C'
binding_indices=metadata_age[metadata_age['subfam']==subfamily_of_interest].index
znf808_sub = znf808[binding_indices]
alignemnt_sub = alignment_filtered[binding_indices]
cov_znf808_mer11c=mapper.normalise(alignment=alignemnt_sub, mapped_data=znf808_sub, method='perc_coverage')
#%%
test_set=metadata_age_sorted.reset_index()
#%%
test_set[(test_set['subfam']=='MER11A')].index
#%%
test_set[(test_set['subfam']=='MER11B')].index
#%%
test_set[(test_set['subfam']=='MER11C')].index
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [znf808_sorted], 
    alignment=alignment_sorted,
    show_alignment=True, 
    heatmap_color=['Greens','Blues'],
    heatmap_mode='overlay',
    heatmap_title=['ZNF808 ChIP-exo signal peaks on MER11 MSA'],
    heatmap_title_fs = 10, 
    anno_ylabel = 'sequences',
    anno_ylabel_fs=10,
    vlim = [[0,7]], 
    opacity = 0.5, 
    hm_transparency_mode = 'gradient',
    aggregated=True, 
    aggregated_data=[cov_znf808,[cov_znf808_mer11a,cov_znf808_mer11b,cov_znf808_mer11c]], 
    agg_colset=['green',['red','blue','yellow']],
    agg_ylim=[[0,110],[0,110]],
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
        annotation=True, 
    anno_col = ['Blues',['red', 'yellow', 'blue']], 
    #anno_title=['TEA-TIME','subfamily'],
    annotation_data=[age_anno_sorted,subfam_anno_sorted],
    anno_cbar_label=[age_subgroups,['MER11A','MER11B','MER11C']],
    anno_cbar_title=['TEA-TIME','subfamily'],  
    anno_cbar_even_pos=-0.19,
        heatmap_xhighlight=[[940,940],[1475,1475]],
    heatmap_xhighlight_col= ['black','black','black'],
    heatmap_xhighlight_alpha=[0.5,0.5,0.5],

    )
#%%
np.nanmean(cov_znf808_mer11c[600:700])
#%%
bam_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/KZFP-bam_hg38/znf808.sorted.bam'
bam_forward=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_forward', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
bam_forward=bam_forward[noNA_indices]
bam_reverse=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_reverse', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
bam_reverse=bam_reverse[noNA_indices]
mean_forward=mapper.normalise(alignment=alignment_filtered, mapped_data=bam_forward)
mean_reverse=mapper.normalise(alignment=alignment_filtered, mapped_data=bam_reverse)
#%%
mean_min=np.minimum(mean_forward,mean_reverse)
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    heatmap=False,
    show_alignment=False,
    alignment=alignment_filtered,
    #logos=True, 
    aggregated=True, 
    aggregated_data=[[mean_forward,mean_reverse], mean_min], 
    agg_colset=[['red','blue'],['green']],
    ag_opacity=0.5,
    agg_ylim=[[0,1.1],[0,1.1]],
    agg_titles=['forward/\nreverse reads','minimum reads'], 
    agg_titles_fs=10,
    agg_titles_pos=[1.2,0.5],
    agg_plot_title=['ZNF808 ChIP-exo signals on MER11 MSA',None],
    agg_ylabel_ypos=[0.05,None],
    agg_plot_title_fs= 12, 
    agg_ylabel=['signal coverage',None],
    agg_ylabel_fs=10,
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    agg_yhighlight=[[450,800],[600,625]],
    agg_yhighlight_col= ['grey','grey'],
    agg_yhighlight_alpha=[0.1,0.1],
    agg_major_tick=100,
    agg_h=20,
    figsize= [80,40],
    #xlim = [212,232]
    )
#%%
metadata_age_sorted=metadata_age.sort_values(['subfam','znf808_group','te_age'])
sorted_indices= metadata_age_sorted.index
znf808_sorted=znf808[sorted_indices]
znf525_sorted=znf525[sorted_indices]
alignment_sorted = alignment_filtered[sorted_indices]
age_anno_sorted=age_anno[sorted_indices]
znf808_subgroup_anno_sorted = znf808_subgroup_anno[sorted_indices]
znf525_subgroup_anno_sorted = znf525_subgroup_anno[sorted_indices]
subfam_anno_sorted = subfam_anno[sorted_indices]
subgroup_anno_sorted=subgroup_anno[sorted_indices]
phylop_sorted = phylop_447[sorted_indices]

#%%
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [znf808_sorted], 
    alignment=alignment_sorted,
    show_alignment=True, 
    heatmap_color=['Greens'],
    heatmap_mode='overlay', 
    vlim = [[0,7]], 
    opacity = 0.5, 
    heatmap_title=['ZNF808 ChIP-exo signal peaks on THE1C, grouped by ZNF808 signals'],
    heatmap_title_fs = 10, 
       heatmap_xhighlight=[[940,940],[1475,1475]],
    heatmap_xhighlight_col= ['black','black','black'],
    heatmap_xhighlight_alpha=[0.5,0.5,0.5],
    heatmap_xlabel = 'position (bp)',
    heatmap_xlabel_fs=10,
    heatmap_ylabel= 'sequences',
    heatmap_ylabel_fs=10,
        agg_major_tick=100,
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
    anno_col = ['Blues',['red', 'yellow', 'blue']], 
    #anno_title=['TEA-TIME','signals','subfamily'],
    annotation_data=[age_anno_sorted,subfam_anno_sorted],
    anno_cbar_label=[age_subgroups,['MER11A','MER11B','MER11C']],
    anno_cbar_title=['TEA-TIME','signals','subfamily'],  
    anno_cbar_even_pos=-0.19,
    colorbar=True,
    colorbar_steps=[0.1],
    heatmap_title=['phyloP of THE1C, grouped by ZNF808 signals'],
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
binding_indices=metadata_age_filtered[(metadata_age_filtered['znf808_group']=='A')].index
nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['znf808_group']=='B')].index
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
    agg_plot_title=[f'phyloP of {subfamily_of_interest}, grouped by ZNF808 signals',None],
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
    binding_indices = metadata_age_filtered[(metadata_age_filtered['znf808_group']=='A')].index
    nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['znf808_group']=='B')].index
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
age_of_interest=20.19
subfamily_of_interest = 'MER11A'
metadata_age_filtered_=metadata_age[metadata_age['subfam']==subfamily_of_interest]
metadata_age_filtered=metadata_age_filtered_[metadata_age_filtered_['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['znf808_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['znf808_group']=='B')].index
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
    agg_ylim=[[0,5],[-5,0]],
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
#%%
#%%
subfamily_of_interest = 'MER11B'
metadata_age_filtered=metadata_age[metadata_age['subfam']==subfamily_of_interest]
binding_indices=metadata_age_filtered[(metadata_age_filtered['znf808_group']=='A')].index
nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['znf808_group']=='B')].index
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
    agg_plot_title=[f'phyloP of {subfamily_of_interest}, grouped by ZNF808 signals',None],
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
    binding_indices = metadata_age_filtered[(metadata_age_filtered['znf808_group']=='A')].index
    nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['znf808_group']=='B')].index
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
age_of_interest=20.19
subfamily_of_interest = 'MER11B'
metadata_age_filtered_=metadata_age[metadata_age['subfam']==subfamily_of_interest]
metadata_age_filtered=metadata_age_filtered_[metadata_age_filtered_['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['znf808_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['znf808_group']=='B')].index
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
#%%
subfamily_of_interest = 'MER11C'
metadata_age_filtered=metadata_age[metadata_age['subfam']==subfamily_of_interest]
binding_indices=metadata_age_filtered[(metadata_age_filtered['znf808_group']=='A')].index
nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['znf808_group']=='B')].index
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
    agg_ylim=[[0,40],[-40,0]],
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
    agg_plot_title=[f'phyloP of {subfamily_of_interest}, grouped by ZNF808 signals',None],
    agg_plot_title_fs= 12, 
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    agg_plottext=[None, f'TE w/ signals: {len(binding_indices)}\nTE w/o signals: {len(nonbinding_indices)}'],
    agg_plottext_fs=10,
    agg_plottext_pos=[0.99,0.01],
    agg_major_tick=100,
    colorbar=False,
    agg_h=15,
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
    binding_indices = metadata_age_filtered[(metadata_age_filtered['znf808_group']=='A')].index
    nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['znf808_group']=='B')].index
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
age_of_interest=20.19
subfamily_of_interest = 'MER11C'
metadata_age_filtered_=metadata_age[metadata_age['subfam']==subfamily_of_interest]
metadata_age_filtered=metadata_age_filtered_[metadata_age_filtered_['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['znf808_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['znf808_group']=='B')].index
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
    agg_h=15,
    figsize= [80,40],
    #xlim=[304,313],
    #gg_major_tick=1,
    )
#%%
age_of_interest=15.76
subfamily_of_interest = 'MER11C'
metadata_age_filtered_=metadata_age[metadata_age['subfam']==subfamily_of_interest]
metadata_age_filtered=metadata_age_filtered_[metadata_age_filtered_['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['znf808_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['znf808_group']=='B')].index
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
#%%
age_of_interest=9.06
subfamily_of_interest = 'MER11C'
metadata_age_filtered_=metadata_age[metadata_age['subfam']==subfamily_of_interest]
metadata_age_filtered=metadata_age_filtered_[metadata_age_filtered_['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['znf808_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['znf808_group']=='B')].index
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
