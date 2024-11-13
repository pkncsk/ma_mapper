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
import importlib
importlib.reload(mapper)
bed_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/kzfp_peak_bed/hg38_kzfps_combined.bed'
kzfp_df=pd.read_csv(bed_filepath, sep='\t', header=None)
kzfp_df.columns=['chrom','start','end','name','score','strand']
bed_file=kzfp_df[kzfp_df['name'].str.contains('ZNF267')]
znf267=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=False, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
znf267=znf267[noNA_indices]
#%%
import importlib
importlib.reload(mapper)
mean_phylop=mapper.normalise(alignment=alignment_filtered, mapped_data=phylop_447)
cov_znf267=mapper.normalise(alignment=alignment_filtered, mapped_data=znf267, method='perc_coverage')
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [znf267,], 
    alignment=alignment_filtered,
    show_alignment=True, 
    heatmap_color=['Purples'],
    heatmap_mode='overlay', 
    vlim = [[0,7]], 
    hm_opacity = 1, 
    hm_transparency_mode = 'gradient',
    #aggregated=True, 
    #aggregated_data=[cov_znf267], 
    #agg_colset=['purple',],
    #agg_ylim=[[0,70]],
    #agg_titles=['znf267 motif'], 
    #agg_ylabel=['perc_coverage'],
    #annotation=True, 
    #anno_col = ['Blues'], 
    #annotation_data=[age_anno],
    #anno_cbar_label=[age_subgroups],
    #anno_cbar_title=['TEA-TIME'], 
    colorbar=False,
    hm_plot_title = 'ZNF267 ChIP-exo signal peaks on THE1C MSA',
    hm_xlabel = 'position (bp)',
    hm_ylabel = 'sequences',
    )
#%%
bam_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/KZFP-bam_hg38/znf267.sorted.bam'
bam_forward=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_forward', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
bam_forward=bam_forward[noNA_indices]
bam_reverse=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_reverse', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
bam_reverse=bam_reverse[noNA_indices]
mean_forward=mapper.normalise(alignment=alignment_filtered, mapped_data=bam_forward)
mean_reverse=mapper.normalise(alignment=alignment_filtered, mapped_data=bam_reverse)
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [bam_forward,bam_reverse], 
    alignment=alignment_filtered,
    #show_alignment=True, 
    heatmap_color=['Blues','Reds'],
    heatmap_mode='overlay', 
    vlim = [[0,0.01]], 
    hm_opacity = 0.9, 
    hm_transparency_mode = 'gradient', 
    aggregated=True, 
    aggregated_data=[[mean_forward,mean_reverse]], 
    agg_colset=[['red','blue']],
    ag_opacity=0.5,
    agg_ylim=[[0,0.1]],
    agg_titles=['znf267 ChIP-exo','znf267 ChIP-exo reverse'], 
    agg_ylabel=['signal coverage','signal coverage'],
    agg_xlabel = 'position (bp)',
    #annotation=True, 
    hm_plot_title = 'ZNF267 ChIP-exo signals on THE1C MSA',
    hm_xlabel = 'position (bp)',
    hm_ylabel = 'sequences',
    #anno_col = ['Blues'], 
    #annotation_data=[age_anno],
    #anno_cbar_label=[age_subgroups],
    #anno_cbar_title=['TEA-TIME'], 
    #colorbar=False,
    )
# %%
#find peaks
import scipy
peaks, _ = scipy.signal.find_peaks(cov_znf267, width = 6)
#%%
#right,left
#highest_peak_index = peaks[np.argmax(mean_znf808[peaks])]
binding_indices = np.unique(np.where(znf267[:, peaks] != 0)[0])
nonbinding_indices=list(set(np.arange(znf267.shape[0])) - set(binding_indices))
#%%
#%%
metadata_age['motif_group'] = 'No Group'
metadata_age.loc[metadata_age.index.isin(binding_indices), 'motif_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices), 'motif_group'] = 'B'


#%%
metadata_sorted=metadata_age.sort_values(['motif_group','te_age'])
sorted_indices=metadata_sorted.index
subgroups = np.unique(metadata_sorted['motif_group'].astype(str))
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno=metadata_sorted['motif_group'].map(numerical_subgroup)
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno=metadata_age['motif_group'].map(numerical_subgroup)
subgroup_anno_sorted=subgroup_anno[sorted_indices]
#%%
phylop_sorted=phylop_447[sorted_indices]
znf267_sorted=znf267[sorted_indices]
alignment_sorted=alignment_filtered[sorted_indices]
age_anno_sorted=age_anno[sorted_indices]
#%%
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [znf267_sorted], 
    alignment=alignment_sorted,
    show_alignment=True, 
    heatmap_color=['Purples'],
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
    anno_col = ['Blues',['purple','white']], 
    annotation_data=[age_anno_sorted,subgroup_anno_sorted],
    anno_cbar_label=[age_subgroups,['znf267','no signal']],
    anno_cbar_title=['TEA-TIME','binding signal'], 
    colorbar=True,
    colorbar_steps=[0.1]
    )
#%%
binding_indices=metadata_age[(metadata_age['motif_group']=='A')].index
nonbinding_indices=metadata_age[(metadata_age['motif_group']=='B')].index
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
    agg_ylim=[[-1.5,1.5],[0,5],[-5,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 43.2
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices=metadata_age_filtered[(metadata_age_filtered['motif_group']=='A')].index
nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['motif_group']=='B')].index
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
    aggregated_data=[-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['blue','red'],
    agg_ylim=[[0,5],[-5,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=['-log10p_adj','log10p_adj'],
    colorbar=False,
    #xlim=[304,313],
    #gg_major_tick=1,
    )
#%%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import gridspec
motif = 'AP-1'
#pos_set = np.arange(61, 69)
#pos_set = np.append(pos_set,[176,250,269])
#pos_set = np.arange(224, 232)
#pos_set = np.append(pos_set,[13,51,212,216,363])
#pos_set = np.arange(291, 299)

pos_set = [ 29,  31,  51,  77, 152, 170, 183, 210, 211, 212, 217, 234, 323,327]
fig = plt.figure(figsize=(3, 1 * len(pos_set)))
gs = gridspec.GridSpec(len(pos_set), 1, height_ratios=[1] * len(pos_set), hspace=0)  # Set hspace to 0
bin_edges = np.concatenate(([-np.inf], np.arange(-3, 3.1, 0.1), [np.inf]))  
for idx, pos in enumerate(pos_set):
    ax = fig.add_subplot(gs[idx])  # Create subplot using GridSpec
    ax.hist(phylop_nonbind[:, pos], bins=bin_edges, color='red', alpha=0.1, label=f'w/o {motif}', density=True)
    ax.hist(phylop_bind[:, pos], bins=bin_edges, color='blue', alpha=0.1, label=f'w/ {motif}', density=True)
    sns.kdeplot(phylop_nonbind[:, pos], color='red', ax=ax, label=f'w/o {motif} KDE', linewidth=1)
    sns.kdeplot(phylop_bind[:, pos], color='blue', ax=ax, label=f'w/ {motif} KDE', linewidth=1)
    ax.set_title(f'{pos}', x=0.001, y=1.0, pad=-14, loc='left')
    ax.set_xlim(-3, 3)
    ax.title.set_position([0.1, 0.95])
    #ax.legend(fontsize=5)
ax.set_xlabel('phyloP')  # Set x-label for the last subplot
plt.tight_layout()
#%%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import gridspec
motif = 'AP-1'
#pos_set = np.arange(61, 69)
#pos_set = np.append(pos_set,[176,250,269])
#pos_set = np.arange(224, 232)
#pos_set = np.append(pos_set,[13,51,212,216,363])
#pos_set = np.arange(291, 299)


fig, ax = plt.subplots(figsize=(6, 2))
bin_edges = np.concatenate(([-np.inf], np.arange(-3, 3.1, 0.1), [np.inf]))  
 # Create subplot using GridSpec
ax.hist(np.ravel(phylop_nonbind), bins=bin_edges, color='red', alpha=0.1, label=f'w/o {motif}', density=True)
ax.hist(np.ravel(phylop_bind), bins=bin_edges, color='blue', alpha=0.1, label=f'w/ {motif}', density=True)
sns.kdeplot(np.ravel(phylop_nonbind), color='red', ax=ax, label=f'w/o {motif} KDE', linewidth=1)
sns.kdeplot(np.ravel(phylop_bind), color='blue', ax=ax, label=f'w/ {motif} KDE', linewidth=1)
ax.set_xlim(-3, 3)
#ax.title.set_position([0.1, 0.95])
    #ax.legend(fontsize=5)
ax.set_xlabel('phyloP')  # Set x-label for the last subplot
plt.tight_layout()
#%%
age_of_interest = 29.44
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices=metadata_age_filtered[(metadata_age_filtered['motif_group']=='A')].index
nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['motif_group']=='B')].index
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
    aggregated_data=[-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['blue','red'],
    agg_ylim=[[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=['-log10p_adj','log10p_adj'],
    colorbar=False,
    agg_h=10,
    figsize=[50,20] 
    #xlim=[304,313],
    #gg_major_tick=1,
    )
#%%
age_of_interest = 20.19
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices=metadata_age_filtered[(metadata_age_filtered['motif_group']=='A')].index
nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['motif_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_ylim=[[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=['-log10p_adj','log10p_adj'],
    colorbar=False,
    agg_h=10,
    figsize=[50,20] 
    #xlim=[304,313],
    #gg_major_tick=1,
    )
#%%
age_of_interest = 15.76
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices=metadata_age_filtered[(metadata_age_filtered['motif_group']=='A')].index
nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['motif_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_ylim=[[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=['-log10p_adj','log10p_adj'],
    colorbar=False,
    agg_h=10,
    figsize=[50,20] 
    #xlim=[304,313],
    #gg_major_tick=1,
    )
#%%
age_of_interest = 9.06
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices=metadata_age_filtered[(metadata_age_filtered['motif_group']=='A')].index
nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['motif_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_ylim=[[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=['-log10p_adj','log10p_adj'],
    colorbar=False,
    agg_h=10,
    figsize=[50,20] 
    #xlim=[304,313],
    #gg_major_tick=1,
    )

#%%
age_of_interest = 0
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices=metadata_age_filtered[(metadata_age_filtered['motif_group']=='A')].index
nonbinding_indices=metadata_age_filtered[(metadata_age_filtered['motif_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_ylim=[[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=['-log10p_adj','log10p_adj'],
    colorbar=False,
    agg_h=10,
    figsize=[50,20] 
    #xlim=[304,313],
    #gg_major_tick=1,
    )
#%%
metadata_age_filtered=metadata_age
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='A')&(metadata_age_filtered['left_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 43.2
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='A')&(metadata_age_filtered['left_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 29.44
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='A')&(metadata_age_filtered['left_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 20.19
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='A')&(metadata_age_filtered['left_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 15.76
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='A')&(metadata_age_filtered['left_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 6.7
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='A')&(metadata_age_filtered['left_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 0
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='A')&(metadata_age_filtered['left_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
metadata_age_filtered=metadata_age
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='A')&(metadata_age_filtered['left_group']=='B')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,10],[-15,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 43.2
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='A')&(metadata_age_filtered['left_group']=='B')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
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
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,7.5],[-15,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 29.44
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='A')&(metadata_age_filtered['left_group']=='B')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 20.19
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='A')&(metadata_age_filtered['left_group']=='B')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 15.76
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='A')&(metadata_age_filtered['left_group']=='B')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 9.06
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='A')&(metadata_age_filtered['left_group']=='B')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 6.7
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='A')&(metadata_age_filtered['left_group']=='B')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 0
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='A')&(metadata_age_filtered['left_group']=='B')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
metadata_age_filtered=metadata_age
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-5,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 43.2
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-5,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 29.44
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 20.19
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 15.76
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 9.06
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 6.7
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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
age_of_interest = 0
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
alignemnt_bind=alignment_filtered[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
#stat_v_greater, p_value_greater = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_greater))
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
#stat_v_less, p_value_less = scipy.stats.ttest_ind(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(np.nan_to_num(p_value_less))
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
    agg_colset=['grey','purple','red'],
    agg_ylim=[[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[100,350]],
    agg_yhighlight_col= ['purple'],
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