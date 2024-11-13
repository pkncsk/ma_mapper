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
vcf_dir = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/vcf-gnomad'

vcf=mapper.map_and_overlay(alignment_file, coord_file, vcf_dir, data_format='vcf', vcf_format='gnomad', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
# %%
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/AP-1(bZIP).bed'
ap1=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
#%%
import importlib
importlib.reload(mapper)
mean_vcf=mapper.normalise(alignment=alignment_filtered, mapped_data=vcf)
cov_ap1=mapper.normalise(alignment=alignment_filtered, mapped_data=ap1, method='perc_coverage')
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [ap1,], 
    alignment=alignment_filtered,
    show_alignment=True, 
    heatmap_color=['Greens'],
    heatmap_mode='overlay', 
    vlim = [[0,7]], 
    opacity = 0.5, 
    aggregated=False, 
    aggregated_data=[cov_ap1], 
    agg_colset=['green',],
    agg_ylim=[[0,70]],
    agg_titles=['AP-1 motif'], 
    agg_ylabel=['perc_coverage'],
    colorbar=False,
    )
# %%
#find peaks
import scipy
peaks, _ = scipy.signal.find_peaks(cov_ap1, width = 6)
#highest_peak_index = peaks[np.argmax(mean_znf808[peaks])]
binding_indices = np.unique(np.where(ap1[:, peaks] != 0)[0])
nonbinding_indices=list(set(np.arange(ap1.shape[0])) - set(binding_indices))
binding_indices_306 = np.unique(np.where(ap1[:, 306] != 0)[0])
binding_indices_100 = np.unique(np.where(ap1[:, 100] != 0)[0])
union_100_306=list(set(binding_indices_100)|set(binding_indices_306))
intersect_100_306=list(set(binding_indices_100)&set(binding_indices_306))
outset_100_306 = list(set(binding_indices)-set(union_100_306))
only_100=list(set(binding_indices_100)-set(intersect_100_306))
only_306=list(set(binding_indices_306)-set(intersect_100_306))
sorted_indices= np.concatenate((intersect_100_306,only_306,only_100,outset_100_306,nonbinding_indices))
ap1_sorted=ap1[sorted_indices]
alignment_sorted = alignment_filtered[sorted_indices]
#cov_ap1[98:108]
#array([ 5.16566958, 10.        , 10.04602992, 10.04487401, 10.02411299,
#       10.02641553, 10.02763067, 10.03460208, 10.0196691 ,  4.75575388])
#cov_ap1[304:314]
#array([32.53400143, 65.06432749, 65.30588373, 65.23914814, 65.02491598,
#       65.12005568, 65.09335498, 64.82903001, 64.75171532, 32.51719716])
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [ap1_sorted,], 
    alignment=alignment_sorted,
    show_alignment=True, 
    heatmap_color=['Greens'],
    heatmap_mode='overlay', 
    vlim = [[0,7]], 
    opacity = 0.5, 

    )
#%%

# %%
vcf_bind = vcf[binding_indices]
vcf_nonbind = vcf[nonbinding_indices]
vcf_sorted = vcf[sorted_indices]

#%%
metadata_age['motif_group'] = 'No Group'
metadata_age.loc[metadata_age.index.isin(intersect_100_306), 'motif_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(only_100), 'motif_group'] = 'B'
metadata_age.loc[metadata_age.index.isin(only_306), 'motif_group'] = 'C'
metadata_age.loc[metadata_age.index.isin(outset_100_306), 'motif_group'] = 'D'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices), 'motif_group'] = 'E'
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
    data = [-np.log10(vcf_sorted),], 
    #alignment=alignment_extended, 
    heatmap_color=['Blues'],
    heatmap_mode='overlay', 
    vlim = [[0,0.00001]],  
    opacity = 0.5, 
    annotation=True, 
    anno_col = [['purple','red','blue','yellow','white']], 
    annotation_data=[subgroup_anno_sorted],
    anno_cbar_label=[['100+306','100 only','306 only', 'others', 'no motif']],
    anno_cbar_title=['motif group'], 
    colorbar=True,
    colorbar_steps = [0.1,0.0001], 
    )
# %%
vcf_bind = vcf[binding_indices]
vcf_nonbind = vcf[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
vcf_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
vcf_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(vcf_bind,vcf_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(vcf_bind,vcf_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(vcf_bind,vcf_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(vcf_bind,vcf_nonbind, axis =0,nan_policy='omit', alternative = 'less')
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
    aggregated_data=[-np.log10(mean_vcf),-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[0,5],[0,25],[-3,0]],
    agg_xhighlight=[[98,108],[304,314]],
    agg_xhighlight_col= ['green','green'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_vcf','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    #xlim=[304,313],
    #gg_major_tick=1,
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
    aggregated_data=[mean_vcf,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    #agg_highlight=[[98,108],[304,314]],
    #agg_highlight_col= ['green','green'],
    #agg_highlight_alpha=[0.2,0.2],
    agg_titles=['mean_vcf','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=['-log10_mean','-log10p_adj','log10p_adj'],
    colorbar=False,
    xlim=[304,313],
    agg_major_tick=1,
    )
# %%
vcf_bind = vcf[intersect_100_306]
vcf_nonbind = vcf[nonbinding_indices]

alignemnt_bind=alignment_filtered[intersect_100_306]
vcf_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
vcf_nonbind[alignemnt_nonbind == 0] = np.nan
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(vcf_bind,vcf_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(vcf_bind,vcf_nonbind, axis =0,nan_policy='omit', alternative = 'less')
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
    aggregated_data=[mean_vcf,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,10],[-15,0]],
    agg_yhighlight=[[98,108],[304,314]],
    agg_yhighlight_col= ['green','green'],
    agg_yhighlight_alpha=[0.2,0.2],#agg_h = 5,
    agg_titles=['mean_vcf','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    )
# %%
vcf_bind = vcf[only_100]
vcf_nonbind = vcf[nonbinding_indices]

alignemnt_bind=alignment_filtered[only_100]
vcf_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
vcf_nonbind[alignemnt_nonbind == 0] = np.nan
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(vcf_bind,vcf_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(vcf_bind,vcf_nonbind, axis =0,nan_policy='omit', alternative = 'less')
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
    aggregated_data=[mean_vcf,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-10,0]],
    agg_yhighlight=[[98,108],[304,314]],
    agg_yhighlight_col= ['green','green'],
    agg_yhighlight_alpha=[0.2,0.2],#agg_h = 5,
    agg_titles=['mean_vcf','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    )
# %%
vcf_bind = vcf[only_306]
vcf_nonbind = vcf[nonbinding_indices]

alignemnt_bind=alignment_filtered[only_306]
vcf_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
vcf_nonbind[alignemnt_nonbind == 0] = np.nan
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(vcf_bind,vcf_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(vcf_bind,vcf_nonbind, axis =0,nan_policy='omit', alternative = 'less')
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
    aggregated_data=[mean_vcf,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,25],[-3,0]],
    agg_yhighlight=[[98,108],[304,314]],
    agg_yhighlight_col= ['green','green'],
    agg_yhighlight_alpha=[0.2,0.2],#agg_h = 5,
    agg_titles=['mean_vcf','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    )
# %%
vcf_bind = vcf[outset_100_306]
vcf_nonbind = vcf[nonbinding_indices]

alignemnt_bind=alignment_filtered[outset_100_306]
vcf_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
vcf_nonbind[alignemnt_nonbind == 0] = np.nan
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(vcf_bind,vcf_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(vcf_bind,vcf_nonbind, axis =0,nan_policy='omit', alternative = 'less')
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
    aggregated_data=[mean_vcf,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    #agg_h = 5,
    agg_titles=['mean_vcf','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    )
#%%
binding_indices_306
# %%
vcf_bind = vcf[binding_indices_306]
vcf_nonbind = vcf[nonbinding_indices]

alignemnt_bind=alignment_filtered[binding_indices_306]
vcf_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
vcf_nonbind[alignemnt_nonbind == 0] = np.nan
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(vcf_bind,vcf_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(vcf_bind,vcf_nonbind, axis =0,nan_policy='omit', alternative = 'less')
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
    aggregated_data=[mean_vcf,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,25],[-3,0]],
    #agg_h = 5,
    agg_titles=['mean_vcf','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    logos=True,
    xlim = [304,313],
    agg_major_tick=1,

    )
#%%

vcf_all_306=vcf[binding_indices_306]
vcf_only_306=vcf[only_306]
vcf_only_100=vcf[only_100]
vcf_100_106=vcf[intersect_100_306]
vcf_others = vcf[outset_100_306]
vcf_nonbind = vcf[nonbinding_indices]
#bindingste 
pos = 313
motif306_all_306=vcf_all_306[:, pos]
motif306_only_306=vcf_only_306[:, pos]
motif306_only_100=vcf_only_100[:, pos]
motif306_100_306=vcf_100_106[:, pos]
motif306_others=vcf_others[:, pos]
motif306_nonbind=vcf_nonbind[:, pos]
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
#ax.set_xlabel('vcf from Zoonomia (241 species)')
#ax.set_ylabel('vcf from Cactus (447 species)')

plt.show()
#%%
motif_100=vcf[:, 98:108]
motif_306=vcf[:,304:314]
mean_motif_100=np.nanmean(motif_100,axis=1)
mean_motif_306=np.nanmean(motif_306,axis=1)
#%%
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6,6))
#ax.hist(motif_100, bins=100, alpha=0.2)
ax.hist(motif_306, bins=100, alpha=0.2)
#ax.set_ylim(-1.2,0.5)
#ax.set_xlim(-1.2,0.5)
#ax.set_xlabel('vcf from Zoonomia (241 species)')
#ax.set_ylabel('vcf from Cactus (447 species)')

plt.show()
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [ap1,], 
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
    aggregated_data=[cov_ap1], 
    agg_colset=['green',],
    agg_ylim=[[0,70]],
    #agg_h = 5,
    agg_titles=['AP-1 motif'], 
    agg_ylabel=['perc_coverage'],
    colorbar=False,
    )
# %% 
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [vcf_sorted,ap1_sorted,], 
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
    aggregated_data=[mean_vcf,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,25],[-3,0]],
    #agg_h = 5,
    agg_titles=['mean_vcf','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=True,
    colorbar_steps = [0.1,0.0001])
#%%
te_age_sorted=metadata_age.iloc[np.concatenate((binding_indices,nonbinding_indices))].te_age.fillna(0)
anno_label=metadata_age.fillna(0).te_age.sort_values().unique()