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
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/NFkB-p65.bed'
nfkb=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
#%%
import importlib
importlib.reload(mapper)
mean_phylop=mapper.normalise(alignment=alignment_filtered, mapped_data=phylop_447)
cov_nfkb=mapper.normalise(alignment=alignment_filtered, mapped_data=nfkb, method='perc_coverage')
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [nfkb,], 
    alignment=alignment_filtered,
    show_alignment=True, 
    heatmap_color=['Blues'],
    heatmap_mode='overlay', 
    vlim = [[0,7]], 
    opacity = 0.5, 
    aggregated=False, 
    aggregated_data=[cov_nfkb], 
    agg_colset=['blue',],
    agg_ylim=[[0,70]],
    agg_titles=['NFkB-p65 motif'], 
    agg_ylabel=['perc_coverage'],
    colorbar=False,
    )
# %%
#find peaks
import scipy
peaks, _ = scipy.signal.find_peaks(cov_nfkb, width = 6)
#highest_peak_index = peaks[np.argmax(mean_znf808[peaks])]
binding_indices = np.unique(np.where(nfkb[:, peaks] != 0)[0])
nonbinding_indices=list(set(np.arange(nfkb.shape[0])) - set(binding_indices))
#%%
binding_indices_180 = np.unique(np.where(nfkb[:, 180] != 0)[0])
binding_indices_115 = np.unique(np.where(nfkb[:, 115] != 0)[0])
union_115_180=list(set(binding_indices_115)|set(binding_indices_180))
intersect_115_180=list(set(binding_indices_115)&set(binding_indices_180))
outset_115_180 = list(set(binding_indices)-set(union_115_180))
only_115=list(set(binding_indices_115)-set(intersect_115_180))
only_180=list(set(binding_indices_180)-set(intersect_115_180))
sorted_indices= np.concatenate((intersect_115_180,only_180,only_115,outset_115_180,nonbinding_indices))
nfkb_sorted=nfkb[sorted_indices]
alignment_sorted = alignment_filtered[sorted_indices]
#cov_nfkb[109:122]
#array([2.98871962, 5.2704899 , 5.33816143, 5.38578621, 5.46115083,
#       5.51474875, 5.55621796, 5.54506843, 3.21680922, 9.85347985,
#       5.4895968 , 5.31457146, 2.95250321])
#cov_nfkb[175:190]
#array([ 7.09976798, 31.04366964, 48.24602432, 49.0115803 , 48.93369071,
#       49.32953602, 48.88507719, 47.40760436, 48.38502123, 48.1260117 ,
#       47.94741079, 46.34054798, 31.57959738, 15.09597523,  7.29379686])
#%%
bed_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/homer_known_motif_hg38/Otx2(Homeobox).bed'
tf=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
tf_sorted= tf[sorted_indices]
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [nfkb_sorted,tf_sorted], 
    alignment=alignment_sorted,
    show_alignment=True, 
    heatmap_color=['Blues','Reds'],
    heatmap_mode='overlay', 
    vlim = [[0,7]], 
    opacity = 0.5, 

    )
#%%

# %%
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
phylop_sorted = phylop_447[sorted_indices]

#%%
metadata_age['motif_group'] = 'No Group'
metadata_age.loc[metadata_age.index.isin(intersect_115_180), 'motif_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(only_115), 'motif_group'] = 'B'
metadata_age.loc[metadata_age.index.isin(only_180), 'motif_group'] = 'C'
metadata_age.loc[metadata_age.index.isin(outset_115_180), 'motif_group'] = 'D'
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
plots.plot_experimental(
    data = [phylop_sorted,], 
    alignment=alignment_sorted,
    #show_alignment=True, 
    heatmap_color=[custom_cmap.vlag_r_mpl],
    heatmap_mode='overlay', 
    vlim = [[-0.5,0.5]], 
    opacity = 0.5, 
    annotation=True, 
    anno_col = [['purple','red','blue','yellow','white']], 
    annotation_data=[subgroup_anno_sorted],
    anno_cbar_label=[['115+180','115 only','180 only', 'others', 'no motif']],
    anno_cbar_title=['motif group'], 
    colorbar=True,
    colorbar_steps=[0.1]
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
    alignment=alignment_sorted,
    #=True, 
    aggregated=True, 
    aggregated_data=[mean_phylop,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,10],[-10,0]],
    agg_yhighlight=[[109,122],[175,190]],
    agg_yhighlight_col= ['blue','blue'],
    agg_yhighlight_alpha=[0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    )
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    heatmap=False,
    show_alignment=False,
    alignment=alignemnt_bind,
    logos=True, 
    aggregated=True, 
    aggregated_data=[mean_phylop,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,10],[-10,0]],
    #agg_highlight=[[109,122],[175,190]],
    #agg_highlight_col= ['blue','blue'],
    #agg_highlight_alpha=[0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    xlim=[175,190],
    agg_major_tick=1,
    )
# %%
phylop_bind = phylop_447[intersect_115_180]
phylop_nonbind = phylop_447[nonbinding_indices]

alignemnt_bind=alignment_filtered[intersect_115_180]
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
    agg_yhighlight=[[109,122],[175,190]],
    agg_yhighlight_col= ['blue','blue'],
    agg_yhighlight_alpha=[0.2,0.2],#agg_h = 5,
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    )
# %%
phylop_bind = phylop_447[only_115]
phylop_nonbind = phylop_447[nonbinding_indices]

alignemnt_bind=alignment_filtered[only_115]
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
    agg_ylim=[[-1.5,1.5],[0,3],[-5,0]],
    agg_yhighlight=[[109,122],[175,190]],
    agg_yhighlight_col= ['blue','blue'],
    agg_yhighlight_alpha=[0.2,0.2],#agg_h = 5,
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    )
# %%
phylop_bind = phylop_447[only_180]
phylop_nonbind = phylop_447[nonbinding_indices]

alignemnt_bind=alignment_filtered[only_180]
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
    agg_yhighlight=[[109,122],[175,190]],
    agg_yhighlight_col= ['blue','blue'],
    agg_yhighlight_alpha=[0.2,0.2],#agg_h = 5,
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    )
# %%
phylop_bind = phylop_447[outset_115_180]
phylop_nonbind = phylop_447[nonbinding_indices]

alignemnt_bind=alignment_filtered[outset_115_180]
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
    agg_yhighlight=[[109,122],[175,190]],
    agg_yhighlight_col= ['blue','blue'],
    agg_yhighlight_alpha=[0.2,0.2],#agg_h = 5,
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    )
#%%
binding_indices_180
# %%
phylop_bind = phylop_447[binding_indices_180]
phylop_nonbind = phylop_447[nonbinding_indices]

alignemnt_bind=alignment_filtered[binding_indices_180]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=scipy.stats.false_discovery_control(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=scipy.stats.false_discovery_control(p_value_less)

mean_phylop_bind=mapper.normalise(alignment=alignemnt_bind, mapped_data=phylop_bind)
mean_phylop_nonbind=mapper.normalise(alignment=alignemnt_nonbind, mapped_data=phylop_nonbind)
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
    agg_ylim=[[-1.5,1.5],[0,10],[-15,0]],
    #agg_h = 5,
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    logos=True,
    xlim = [175,190],
    agg_major_tick=1,
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    )
#%%

phyloP_all_180=phylop_447[binding_indices_180]
phyloP_only_180=phylop_447[only_180]
phyloP_only_115=phylop_447[only_115]
phyloP_115_106=phylop_447[intersect_115_180]
phylop_others = phylop_447[outset_115_180]
phylop_nonbind = phylop_447[nonbinding_indices]
#bindingste 
pos = 180
motif180_all_180=phyloP_all_180[:, pos]
motif180_only_180=phyloP_only_180[:, pos]
motif180_only_115=phyloP_only_115[:, pos]
motif180_115_180=phyloP_115_106[:, pos]
motif180_others=phylop_others[:, pos]
motif180_nonbind=phylop_nonbind[:, pos]
#filter nan
motif180_all_180_filtered = motif180_all_180[~np.isnan(motif180_all_180)]
motif180_only_180_filtered = motif180_only_180[~np.isnan(motif180_only_180)]
motif180_only_115_filtered = motif180_only_115[~np.isnan(motif180_only_115)]
motif180_115_180_filtered = motif180_115_180[~np.isnan(motif180_115_180)]
motif180_others_filtered = motif180_others[~np.isnan(motif180_others)]
motif180_nonbind_filtered = motif180_nonbind[~np.isnan(motif180_nonbind)]
#find average
#mean_motif180_all_180 = np.nanmean(motif180_all_180,axis=1)
#mean_motif180_only_180 = np.nanmean(motif180_only_180,axis=1)
#mean_motif180_only_115 = np.nanmean(motif180_only_115,axis=1)
#mean_motif180_115_180 = np.nanmean(motif180_115_180,axis=1)
#mean_motif180_others = np.nanmean(motif180_others,axis=1)
#mean_motif180_nonbind = np.nanmean(motif180_nonbind,axis=1)
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6,6))
ax.violinplot([motif180_all_180_filtered,motif180_only_180_filtered,motif180_only_115_filtered,motif180_115_180_filtered,motif180_others_filtered,motif180_nonbind_filtered])
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
    data = [nfkb,], 
    alignment=alignment_filtered,
    show_alignment=True, 
    heatmap_color=['Blues'],
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
    aggregated_data=[cov_nfkb], 
    agg_colset=['blue',],
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
    data = [phylop_sorted,nfkb_sorted,], 
    alignment=alignment_sorted,
    #show_alignment=True, 
    heatmap_color=[custom_cmap.vlag_r_mpl,'Blues'],
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
    colorbar_steps = [0.1,0.0001],  )
#%%
te_age_sorted=metadata_age.iloc[np.concatenate((binding_indices,nonbinding_indices))].te_age.fillna(0)
anno_label=metadata_age.fillna(0).te_age.sort_values().unique()