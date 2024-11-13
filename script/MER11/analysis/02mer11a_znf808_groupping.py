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
subfamily = 'MER11A'
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

#%%
coord_file=mapper.extract_metadata_from_alignment(alignment_file)
#%%
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/zoonomia447/hg38.phyloP447way.bw'
phylop_447=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
# %%
bed_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/kzfp_peak_bed/hg38_kzfps_combined.bed'
kzfp_df=pd.read_csv(bed_filepath, sep='\t', header=None)
kzfp_df.columns=['chrom','start','end','name','score','strand']
bed_file=kzfp_df[kzfp_df['name'].str.contains('ZNF808')]
znf808=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap = False, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
#%%
import importlib
importlib.reload(mapper)
mean_phylop=mapper.normalise(alignment=alignment_filtered, mapped_data=phylop_447)
cov_znf808=mapper.normalise(alignment=alignment_filtered, mapped_data=znf808, method='perc_coverage')
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [znf808,], 
    alignment=alignment_filtered,
    show_alignment=True, 
    heatmap_color=['Greens'],
    heatmap_mode='overlay', 
    vlim = [[0,7]], 
    opacity = 0.5, 
    aggregated=True, 
    aggregated_data=[cov_znf808], 
    agg_colset=['green',],
    agg_ylim=[[0,100]],
    agg_titles=['ZNF808 motif'], 
    agg_ylabel=['perc_coverage'],
    colorbar=False,
    agg_major_tick=100
    )
# %%
#find peaks
import scipy
peaks, _ = scipy.signal.find_peaks(cov_znf808, width = 6)
#%%
#highest_peak_index = peaks[np.argmax(mean_znf808[peaks])]
binding_indices = np.unique(np.where(znf808[:, peaks] != 0)[0])
nonbinding_indices=list(set(np.arange(znf808.shape[0])) - set(binding_indices))
sorted_indices= np.concatenate((binding_indices,nonbinding_indices))
znf808_sorted=znf808[sorted_indices]
alignment_sorted = alignment_filtered[sorted_indices]
#cov_znf808[98:108]
#array([ 5.16566958, 10.        , 10.04602992, 10.04487401, 10.02411299,
#       10.02641553, 10.02763067, 10.03460208, 10.0196691 ,  4.75575388])
#cov_znf808[304:314]
#array([32.53400143, 65.06432749, 65.30588373, 65.23914814, 65.02491598,
#       65.12005568, 65.09335498, 64.82903001, 64.75171532, 32.51719716])
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [znf808_sorted,], 
    alignment=alignment_sorted,
    show_alignment=True, 
    heatmap_color=['Greens'],
    heatmap_mode='overlay', 
    vlim = [[0,7]], 
    opacity = 0.5, 
    agg_major_tick=100
    )
#%%

# %%
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
phylop_sorted = phylop_447[sorted_indices]

#%%
metadata_age['motif_group'] = 'No Group'
metadata_age.loc[metadata_age.index.isin(binding_indices), 'motif_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices), 'motif_group'] = 'B'
subgroups = np.unique(metadata_age['motif_group'].astype(str))
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno=metadata_age['motif_group'].map(numerical_subgroup)
subgroup_anno_sorted=subgroup_anno.reindex(sorted_indices)
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(heatmap=True,
    data = [phylop_sorted,], 
    #alignment=alignment_extended, 
    heatmap_color=[custom_cmap.vlag_r_mpl,],
    heatmap_mode='overlay', 
    vlim = [[-0.5,0.5],[0,0.0005],],  
    opacity = 0.5, 
    annotation=True, 
    anno_col = [['green','white'],['red', 'yellow', 'blue']], 
    annotation_data=[subgroup_anno_sorted],
    anno_cbar_label=[['w/ znf808 binding','w/o znf808 binding'],['MER11A','MER11B','MER11C']],
    anno_cbar_title=['motif group','subfamily'], 
    colorbar=True,
    colorbar_steps = [0.1,0.0001], 
    agg_major_tick=100
    )
# %%
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
    agg_ylim=[[-1.5,1.5],[0,30],[-15,0]],
    agg_yhighlight=[[450,800]],
    agg_yhighlight_col= ['green','green'],
    agg_yhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    #xlim=[304,313],
    agg_major_tick=100,
    )
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
    agg_ylim=[[-1.5,1.5],[0,25],[-3,0]],
    #agg_highlight=[[98,108],[304,314]],
    #agg_highlight_col= ['green','green'],
    #agg_highlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    xlim=[304,313],
    agg_major_tick=1,
    )
# %%
phylop_bind = phylop_447[intersect_100_306]
phylop_nonbind = phylop_447[nonbinding_indices]

alignemnt_bind=alignment_filtered[intersect_100_306]
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
    agg_yhighlight=[[98,108],[304,314]],
    agg_yhighlight_col= ['green','green'],
    agg_yhighlight_alpha=[0.2,0.2],#agg_h = 5,
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    )
# %%
phylop_bind = phylop_447[only_100]
phylop_nonbind = phylop_447[nonbinding_indices]

alignemnt_bind=alignment_filtered[only_100]
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
    agg_ylim=[[-1.5,1.5],[0,3],[-10,0]],
    agg_highlight=[[98,108],[304,314]],
    agg_highlight_col= ['green','green'],
    agg_highlight_alpha=[0.2,0.2],#agg_h = 5,
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    )
# %%
phylop_bind = phylop_447[only_306]
phylop_nonbind = phylop_447[nonbinding_indices]

alignemnt_bind=alignment_filtered[only_306]
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
    agg_ylim=[[-1.5,1.5],[0,25],[-3,0]],
    agg_highlight=[[98,108],[304,314]],
    agg_highlight_col= ['green','green'],
    agg_highlight_alpha=[0.2,0.2],#agg_h = 5,
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    )
# %%
phylop_bind = phylop_447[outset_100_306]
phylop_nonbind = phylop_447[nonbinding_indices]

alignemnt_bind=alignment_filtered[outset_100_306]
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
    #agg_h = 5,
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    )
#%%
binding_indices_306
# %%
phylop_bind = phylop_447[binding_indices_306]
phylop_nonbind = phylop_447[nonbinding_indices]

alignemnt_bind=alignment_filtered[binding_indices_306]
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
    alignment=alignemnt_bind,
    aggregated=True, 
    aggregated_data=[mean_phylop,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,25],[-3,0]],
    #agg_h = 5,
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    logos=True,
    xlim = [304,313],
    agg_major_tick=1,

    )
#%%

phyloP_all_306=phylop_447[binding_indices_306]
phyloP_only_306=phylop_447[only_306]
phyloP_only_100=phylop_447[only_100]
phyloP_100_106=phylop_447[intersect_100_306]
phylop_others = phylop_447[outset_100_306]
phylop_nonbind = phylop_447[nonbinding_indices]
#bindingste 
pos = 313
motif306_all_306=phyloP_all_306[:, pos]
motif306_only_306=phyloP_only_306[:, pos]
motif306_only_100=phyloP_only_100[:, pos]
motif306_100_306=phyloP_100_106[:, pos]
motif306_others=phylop_others[:, pos]
motif306_nonbind=phylop_nonbind[:, pos]
#filter nan
motif306_all_306_filtered = motif306_all_306[~np.isnan(motif306_all_306)]
motif306_only_306_filtered = motif306_only_306[~np.isnan(motif306_only_306)]
motif306_only_100_filtered = motif306_only_100[~np.isnan(motif306_only_100)]
motif306_100_306_filtered = motif306_100_306[~np.isnan(motif306_100_306)]
motif306_others_filtered = motif306_others[~np.isnan(motif306_others)]
motif306_nonbind_filtered = motif306_nonbind[~np.isnan(motif306_nonbind)]
#find average
#mean_motif306_all_306 = np.nanmean(motif306_all_306,axis=1)
#mean_motif306_only_306 = np.nanmean(motif306_only_306,axis=1)
#mean_motif306_only_100 = np.nanmean(motif306_only_100,axis=1)
#mean_motif306_100_306 = np.nanmean(motif306_100_306,axis=1)
#mean_motif306_others = np.nanmean(motif306_others,axis=1)
#mean_motif306_nonbind = np.nanmean(motif306_nonbind,axis=1)
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6,6))
ax.violinplot([motif306_all_306_filtered,motif306_only_306_filtered,motif306_only_100_filtered,motif306_100_306_filtered,motif306_others_filtered,motif306_nonbind_filtered],showextrema=False)
#ax.set_ylim(-1.2,0.5)
#ax.set_xlim(-1.2,0.5)
#ax.set_xlabel('phyloP from Zoonomia (241 species)')
#ax.set_ylabel('phyloP from Cactus (447 species)')

plt.show()
#%%
motif_100=phylop_447[:, 98:108]
motif_306=phylop_447[:,304:314]
mean_motif_100=np.nanmean(motif_100,axis=1)
mean_motif_306=np.nanmean(motif_306,axis=1)
#%%
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6,6))
#ax.hist(motif_100, bins=100, alpha=0.2)
ax.hist(motif_306, bins=100, alpha=0.2)
#ax.set_ylim(-1.2,0.5)
#ax.set_xlim(-1.2,0.5)
#ax.set_xlabel('phyloP from Zoonomia (241 species)')
#ax.set_ylabel('phyloP from Cactus (447 species)')

plt.show()
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [znf808,], 
    alignment=alignment_filtered,
    show_alignment=True, 
    heatmap_color=['Greens'],
    heatmap_mode='overlay', 
    vlim = [[0,7]], 
    opacity = 0.5, 
    #annotation=True, 
    #anno_col = ['Blues'], 
    #annotation_data=[te_age_sorted],
    #anno_cbar_label=[anno_label],
    #anno_title=['age'],
    #anno_cbar_title=['MYA'], 
    aggregated=False, 
    aggregated_data=[cov_znf808], 
    agg_colset=['green',],
    agg_ylim=[[0,70]],
    #agg_h = 5,
    agg_titles=['ZNF808 motif'], 
    agg_ylabel=['perc_coverage'],
    colorbar=False,
    )
# %% 
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [phylop_sorted,znf808_sorted,], 
    alignment=alignment_sorted,
    #show_alignment=True, 
    heatmap_color=[custom_cmap.vlag_r_mpl,'Greens'],
    heatmap_mode='spread_horizontal', 
    vlim = [[-0.5,0.5],[0,0.0005],[0,0.0005]], 
    opacity = 0.5, 
    #annotation=True, 
    #anno_col = ['Blues'], 
    #annotation_data=[te_age_sorted],
    #anno_cbar_label=[anno_label],
    #anno_title=['age'],
    #anno_cbar_title=['MYA'], 
    aggregated=True, 
    aggregated_data=[mean_phylop,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,25],[-3,0]],
    #agg_h = 5,
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=True,
    colorbar_steps = [0.1,0.0001])
#%%
te_age_sorted=metadata_age.iloc[np.concatenate((binding_indices,nonbinding_indices))].te_age.fillna(0)
anno_label=metadata_age.fillna(0).te_age.sort_values().unique()