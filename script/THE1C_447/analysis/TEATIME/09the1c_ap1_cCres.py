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
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/AP-1(bZIP).bed'
ap1=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
ap1=ap1[noNA_indices]
#%%
#%% there are more than one approach, this one is just convinient with the pipeline
bed_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/phenotype_tracks/encodeCcreCombined.bed'
ccres_df=pd.read_csv(bed_file, sep='\t', header=None)
ccres_df.columns = ['chrom','start','end','name','score','strand','thickStart','thickEnd','reserved','ccre','encodeLabel','zScore','ucscLabel','accessionLabel','description']
ccres=mapper.map_and_overlay(alignment_file, coord_file, ccres_df, data_format='bed', custom_id=True, strand_overlap=False, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
#%%
import importlib
importlib.reload(mapper)
mean_phylop=mapper.normalise(alignment=alignment_filtered, mapped_data=phylop_447)
cov_ap1=mapper.normalise(alignment=alignment_filtered, mapped_data=ap1, method='perc_coverage')
cov_ccres=mapper.normalise(alignment=alignment_filtered, mapped_data=ccres, method='perc_coverage')
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
    aggregated=True, 
    aggregated_data=[cov_ap1], 
    agg_colset=['green',],
    agg_ylim=[[0,70]],
    agg_titles=['ap1 motif'], 
    agg_ylabel=['perc_coverage'],
    annotation=True, 
    anno_col = ['Blues'], 
    annotation_data=[age_anno],
    anno_cbar_label=[age_subgroups],
    anno_cbar_title=['TEA-TIME'], 
    colorbar=False,
    )
# %%
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
metadata_age['right_group'] = 'No Group'
metadata_age.loc[metadata_age.index.isin(binding_indices_right), 'right_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices_right), 'right_group'] = 'B'
metadata_age['left_group'] = 'No Group'
metadata_age.loc[metadata_age.index.isin(binding_indices_left), 'left_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices_left), 'left_group'] = 'B'
#%%
#find peaks
import scipy
peaks, _ = scipy.signal.find_peaks(cov_ccres, width = 6)
#%%
#right,left
#highest_peak_index = peaks[np.argmax(mean_znf808[peaks])]
annot_indices = np.unique(np.where(ccres[:, peaks] != 0)[0])
nonannot_indices=list(set(np.arange(ccres.shape[0])) - set(annot_indices))
#%%
metadata_age['ccres_group'] = 'No Group'
metadata_age.loc[metadata_age.index.isin(annot_indices), 'ccres_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(nonannot_indices), 'ccres_group'] = 'B'

#%%
metadata_sorted=metadata_age.sort_values(['right_group','left_group','te_age','ccres_group'])
sorted_indices=metadata_sorted.index
subgroups = np.unique(metadata_sorted['ccres_group'].astype(str))
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno_ccres=metadata_sorted['ccres_group'].map(numerical_subgroup)
#%%
phylop_sorted=phylop_447[sorted_indices]
ap1_sorted=ap1[sorted_indices]
alignment_sorted=alignment_filtered[sorted_indices]
age_anno_sorted=age_anno[sorted_indices]
ccres_sorted = ccres[sorted_indices]
#%%

#%%
binding_indices=metadata_sorted[(metadata_sorted['right_group']=='A')|(metadata_sorted['left_group']=='A')].index
nonbinding_indices=metadata_sorted[(metadata_sorted['right_group']=='B')&(metadata_sorted['left_group']=='B')].index
only_right = metadata_sorted[(metadata_sorted['right_group']=='A')&(metadata_sorted['left_group']=='B')].index
only_left = metadata_sorted[(metadata_sorted['right_group']=='B')&(metadata_sorted['left_group']=='A')].index

intersect_right_left = metadata_sorted[(metadata_sorted['right_group']=='A')&(metadata_sorted['left_group']=='A')].index

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
subgroup_anno_ccres_sorted = subgroup_anno_ccres[sorted_indices]
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [ap1_sorted], 
    alignment=alignment_sorted,
    show_alignment=True, 
    heatmap_color=['Greens'],
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
    opacity = 1.0, 
    annotation=True, 
    anno_col = ['Blues',['purple','red','blue','white'],['yellow','purple']], 
    annotation_data=[age_anno_sorted,subgroup_anno_sorted,subgroup_anno_ccres_sorted],
    anno_cbar_label=[age_subgroups,['both','right only','left only','no motif/others'],['cCRES','NA']],
    anno_cbar_title=['TEA-TIME','motif group', 'cCRES'], 
    colorbar=True,
    colorbar_steps=[0.1]
    )



# %%
age_of_interest = 43.2
metadata_age_filtered=metadata_age[metadata_age['te_age']==age_of_interest]
binding_indices = metadata_age_filtered[((metadata_age_filtered['right_group']=='A')|(metadata_age_filtered['left_group']=='A'))&(metadata_age_filtered['ccres_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']=='B')&(metadata_age_filtered['left_group']=='B')&(metadata_age_filtered['ccres_group']=='B')].index
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
    agg_ylim=[[-1.5,1.5],[0,20],[-20,0]],
    agg_yhighlight=[[98,108],[304,314]],
    agg_yhighlight_col= ['green','green'],
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
binding_indices = metadata_age_filtered[(metadata_age_filtered['ccres_group']=='A')].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['ccres_group']=='B')].index
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
    alignment=alignment_sorted,
    #logos=True, 
    aggregated=True, 
    aggregated_data=[mean_phylop_bind,mean_phylop_nonbind,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[-1.5,1.5],[0,3],[-3,0]],
    agg_yhighlight=[[98,108],[304,314]],
    agg_yhighlight_col= ['green','green'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP w/ motif','mean_phyloP w/o motif','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    #xlim=[304,313],
    #gg_major_tick=1,
    agg_h=10
    )
# %%
