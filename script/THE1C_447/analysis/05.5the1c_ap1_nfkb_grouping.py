#%%
import sys
from turtle import width
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
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
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/AP-1(bZIP).bed'
ap1=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/NFkB-p65.bed'
nfkb=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
#%%

#%%
import importlib
importlib.reload(mapper)
mean_phylop=mapper.normalise(alignment=alignment_filtered, mapped_data=phylop_447)
cov_ap1=mapper.normalise(alignment=alignment_filtered, mapped_data=ap1, method='perc_coverage')
cov_nfkb=mapper.normalise(alignment=alignment_filtered, mapped_data=nfkb, method='perc_coverage')

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

# %%
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
phylop_sorted = phylop_447[sorted_indices]

#%%
metadata_age['ap1_group'] = 'No Group'
metadata_age.loc[metadata_age.index.isin(intersect_100_306), 'ap1_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(only_100), 'ap1_group'] = 'B'
metadata_age.loc[metadata_age.index.isin(only_306), 'ap1_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(outset_100_306), 'ap1_group'] = 'B'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices), 'ap1_group'] = 'B'
ap1_subgroups = np.unique(metadata_age['ap1_group'].astype(str))
ap1_numerical_subgroup = {subgroup: num for num, subgroup in enumerate(ap1_subgroups)}
ap1_subgroup_anno=metadata_age['ap1_group'].map(ap1_numerical_subgroup)
ap_1_subgroup_anno_sorted=ap1_subgroup_anno.reindex(sorted_indices)
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
# %%
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]
phylop_sorted = phylop_447[sorted_indices]

#%%
metadata_age['nfkb_group'] = 'No Group'
metadata_age.loc[metadata_age.index.isin(intersect_115_180), 'nfkb_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(only_115), 'nfkb_group'] = 'B'
metadata_age.loc[metadata_age.index.isin(only_180), 'nfkb_group'] = 'A'
metadata_age.loc[metadata_age.index.isin(outset_115_180), 'nfkb_group'] = 'B'
metadata_age.loc[metadata_age.index.isin(nonbinding_indices), 'nfkb_group'] = 'B'
nfkb_subgroups = np.unique(metadata_age['nfkb_group'].astype(str))
nfkb_numerical_subgroup = {subgroup: num for num, subgroup in enumerate(nfkb_subgroups)}
nfkb_subgroup_anno=metadata_age['nfkb_group'].map(nfkb_numerical_subgroup)
nfkb_subgroup_anno_sorted=nfkb_subgroup_anno.reindex(sorted_indices)
#%%
metadata_sorted=metadata_age.sort_values(['ap1_group','nfkb_group'])
final_sort_indices=metadata_sorted.index
nfkb_subgroups = np.unique(metadata_age['nfkb_group'].astype(str))
nfkb_numerical_subgroup = {subgroup: num for num, subgroup in enumerate(nfkb_subgroups)}
nfkb_subgroup_anno=metadata_sorted['nfkb_group'].map(nfkb_numerical_subgroup)
ap1_subgroups = np.unique(metadata_age['ap1_group'].astype(str))
ap1_numerical_subgroup = {subgroup: num for num, subgroup in enumerate(ap1_subgroups)}
ap1_subgroup_anno=metadata_sorted['ap1_group'].map(ap1_numerical_subgroup)
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
    vlim = [[0,7],[0,7]], 
    opacity = 0.5, 
    aggregated=True, 
    aggregated_data=[cov_ap1, cov_nfkb], 
    agg_colset=['green','blue'],
    agg_ylim=[[0,70],[0,70]],
    agg_titles=['AP-1 motif','NFkB-p65 motif'], 
    agg_ylabel=['perc_coverage', 'perc_coverage'],
    colorbar=False,
    )
#%%
phylop_sorted=phylop_447[final_sort_indices]
ap1_sorted=ap1[final_sort_indices]
nfkb_sorted=nfkb[final_sort_indices]
alignemnt_sorted=alignment_filtered[final_sort_indices]
tf_sorted=tf[final_sort_indices]
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [ap1_sorted,nfkb_sorted,tf_sorted], 
    alignment=alignment_sorted,
    show_alignment=True, 
    heatmap_color=['Greens','Blues','Reds'],
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
    anno_col = [['green','white'],['blue','white']], 
    annotation_data=[ap1_subgroup_anno,nfkb_subgroup_anno],
    anno_cbar_label=[['AP-1_306','no AP-1-306'],['NFkB_180','no NFkB_180']],
    anno_cbar_title=['AP-1','NFkB'], 
    colorbar=True,
    colorbar_steps=[0.1]
    )
# %%
binding_indices=metadata_sorted[(metadata_sorted['ap1_group']=='A')&(metadata_sorted['nfkb_group']=='A')].index
nonbinding_indices=metadata_sorted[(metadata_sorted['ap1_group']=='B')&(metadata_sorted['nfkb_group']=='B')].index
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
    agg_ylim=[[-1.5,1.5],[0,20],[-5,0]],
    agg_yhighlight=[[175,190],[304,314]],
    agg_yhighlight_col= ['blue','green'],
    agg_yhighlight_alpha=[0.2,0.2,0.2],
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
# %%
binding_indices=metadata_sorted[(metadata_sorted['ap1_group']=='A')&(metadata_sorted['nfkb_group']=='B')].index
nonbinding_indices=metadata_sorted[(metadata_sorted['ap1_group']=='B')&(metadata_sorted['nfkb_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]

alignemnt_bind=alignment_filtered[binding_indices]
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
    agg_ylim=[[-1.5,1.5],[0,15],[-3,0]],

    agg_yhighlight=[[109,122],[175,190],[304,314]],
    agg_yhighlight_col= ['blue','blue','green'],
    agg_yhighlight_alpha=[0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    )
# %%
binding_indices=metadata_sorted[(metadata_sorted['ap1_group']=='B')&(metadata_sorted['nfkb_group']=='A')].index
nonbinding_indices=metadata_sorted[(metadata_sorted['ap1_group']=='B')&(metadata_sorted['nfkb_group']=='B')].index
phylop_bind = phylop_447[binding_indices]
phylop_nonbind = phylop_447[nonbinding_indices]

alignemnt_bind=alignment_filtered[binding_indices]
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
    agg_ylim=[[-1.5,1.5],[0,7.5],[-3,0]],
    agg_yhighlight=[[109,122],[175,190],[304,314]],
    agg_yhighlight_col= ['blue','blue','green'],
    agg_yhighlight_alpha=[0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
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
    alignment=alignemnt_bind,
    aggregated=True, 
    aggregated_data=[mean_phylop,-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,6],[-5,0]],
    #agg_h = 5,
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    logos=True,
    xlim = [40,60],
    agg_major_tick=1,

    )
# %%
alignment_filtered[:,40:60]
# %%
# Mapping dictionary
mapping = {0: '-', 1: 'A', 2: 'C', 3: 'T', 4: 'G'}

# Convert the matrix to FASTA format
def matrix_to_fasta(matrix, mapping):
    fasta_str = ""
    for i, row in enumerate(matrix):
        # Convert each row to its corresponding character sequence
        sequence = ''.join(mapping[int(num)] for num in row)
        # Add the FASTA format header and sequence
        fasta_str += f">Sequence_{i+1}\n{sequence}\n"
    return fasta_str

# Convert the matrix and print the FASTA formatted string
fasta_formatted_str = matrix_to_fasta(alignemnt_bind[:,40:60], mapping)
print(fasta_formatted_str)

# %%
