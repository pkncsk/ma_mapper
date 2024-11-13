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
#%%
atac_table = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/phenotype_tracks/encode_v3_score_matrix_norm_class_kzfps_embed_all_clust_TSS.csv'
atac_df=pd.read_csv(atac_table, sep=',', low_memory=False, index_col = 0)
#%%
#%%
#atac_df_bed = atac_df[['chrom','start','end','name']]
#atac_df_bed['score']=10
#atac_df_bed['strand'] = '.'
atac_df_bed = pd.read_csv('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/phenotype_tracks/atac_coord_hg38.bed', sep='\t', header = None)
atac=mapper.map_and_overlay(alignment_file, coord_file, atac_df_bed, data_format='bed', custom_id=True, strand_overlap=False, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
atac=atac[noNA_indices]
#%%
bed_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/phenotype_tracks/encodeCcreCombined.bed'
ccres_df=pd.read_csv(bed_file, sep='\t', header=None)
ccres_df.columns = ['chrom','start','end','name','score','strand','thickStart','thickEnd','reserved','ccre','encodeLabel','zScore','ucscLabel','accessionLabel','description']
ccres=mapper.map_and_overlay(alignment_file, coord_file, ccres_df, data_format='bed', custom_id=True, strand_overlap=False, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
#%%
#%%
bed_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/phenotype_tracks/ZNF808-active.bed'
znf808active_df=pd.read_csv(bed_file, sep='\t', header=None)

#%%
import importlib
importlib.reload(mapper)
mean_phylop=mapper.normalise(alignment=alignment_filtered, mapped_data=phylop_447)
cov_atac=mapper.normalise(alignment=alignment_filtered, mapped_data=atac, method='perc_coverage')
cov_ap1=mapper.normalise(alignment=alignment_filtered, mapped_data=ap1, method='perc_coverage')
cov_ccres=mapper.normalise(alignment=alignment_filtered, mapped_data=ccres, method='perc_coverage')
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
metadata_age['right_group'] = 1
metadata_age.loc[metadata_age.index.isin(binding_indices_right), 'right_group'] = 0
metadata_age.loc[metadata_age.index.isin(nonbinding_indices_right), 'right_group'] = 1
metadata_age['left_group'] = 1
metadata_age.loc[metadata_age.index.isin(binding_indices_left), 'left_group'] = 0
metadata_age.loc[metadata_age.index.isin(nonbinding_indices_left), 'left_group'] = 1
#%%
#find peaks
import scipy
peaks, _ = scipy.signal.find_peaks(cov_atac, width = 6)
#%%
#right,left
#highest_peak_index = peaks[np.argmax(mean_znf808[peaks])]
annot_indices = np.unique(np.where(atac[:, peaks] != 0)[0])
nonannot_indices=list(set(np.arange(atac.shape[0])) - set(annot_indices))
#%%
metadata_age['atac_group'] = 0
metadata_age.loc[metadata_age.index.isin(annot_indices), 'atac_group'] = 0
metadata_age.loc[metadata_age.index.isin(nonannot_indices), 'atac_group'] = 1

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
metadata_age['ccres_group'] = 0
metadata_age.loc[metadata_age.index.isin(annot_indices), 'ccres_group'] = 0
metadata_age.loc[metadata_age.index.isin(nonannot_indices), 'ccres_group'] = 1

#%%
metadata_sorted=metadata_age.sort_values(['right_group','left_group','te_age','atac_group','ccres_group'])
sorted_indices=metadata_sorted.index
#subgroups = np.unique(metadata_sorted['atac_group'].astype(str))
#numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
#subgroup_anno_atac=metadata_sorted['atac_group'].map(numerical_subgroup)
#subgroups = np.unique(metadata_sorted['ccres_group'].astype(str))
#numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
#subgroup_anno_ccres=metadata_sorted['ccres_group'].map(numerical_subgroup)
subgroup_anno_atac=metadata_sorted['atac_group']
subgroup_anno_ccres=metadata_sorted['ccres_group']
#%%
phylop_sorted=phylop_447[sorted_indices]
ap1_sorted=ap1[sorted_indices]
alignment_sorted=alignment_filtered[sorted_indices]
age_anno_sorted=age_anno[sorted_indices]
atac_sorted = atac[sorted_indices]
ccres_sorted = ccres[sorted_indices]
#%%

#%%
binding_indices=metadata_sorted[(metadata_sorted['right_group']==0)|(metadata_sorted['left_group']==0)].index
nonbinding_indices=metadata_sorted[(metadata_sorted['right_group']==1)&(metadata_sorted['left_group']==1)].index
only_right = metadata_sorted[(metadata_sorted['right_group']==0)&(metadata_sorted['left_group']==1)].index
only_left = metadata_sorted[(metadata_sorted['right_group']==1)&(metadata_sorted['left_group']==0)].index

intersect_right_left = metadata_sorted[(metadata_sorted['right_group']==0)&(metadata_sorted['left_group']==0)].index

#%%
metadata_age['motif_group'] = 0
metadata_age.loc[metadata_age.index.isin(intersect_right_left), 'motif_group'] = 0
metadata_age.loc[metadata_age.index.isin(only_right), 'motif_group'] = 1
metadata_age.loc[metadata_age.index.isin(only_left), 'motif_group'] = 2
metadata_age.loc[metadata_age.index.isin(nonbinding_indices), 'motif_group'] = 3
#subgroups = np.unique(metadata_age['motif_group'].astype(str))
#numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
#subgroup_anno=metadata_age['motif_group'].map(numerical_subgroup)
subgroup_anno=metadata_age['motif_group']
subgroup_anno_sorted=subgroup_anno[sorted_indices]
subgroup_anno_atac_sorted = subgroup_anno_atac[sorted_indices]
subgroup_anno_ccres_sorted = subgroup_anno_ccres[sorted_indices]
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
    vlim = [[-1.0,1.0]], 
    opacity = 1.0, 
    annotation=True, 
    anno_col = ['Blues',['purple','red','blue','white'],['green','white'],['orange','white']], 
    annotation_data=[age_anno_sorted,subgroup_anno_sorted,subgroup_anno_atac_sorted,subgroup_anno_ccres_sorted],
    anno_cbar_label=[age_subgroups,['both','right only','left only','no motif/others'],['ATAC','NA'],['cCRES','NA']],
    anno_cbar_title=['TEA-TIME','motif group', 'ATAC','cCRES'], 
    colorbar=True,
    colorbar_steps=[0.1]
    )


# %%
# %%
metadata_age_filtered=metadata_age
binding_indices = metadata_age_filtered[((metadata_age_filtered['right_group']==0)|(metadata_age_filtered['left_group']==0))&(metadata_age_filtered['atac_group']==0)].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['right_group']==1)&(metadata_age_filtered['left_group']==1)&(metadata_age_filtered['atac_group']==1)].index
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
# %%
#atac_df_bed = atac_df[['chrom','start','end','name']]
#atac_df_bed['score']=10
#atac_df_bed['strand'] = '.'
import pybedtools
atac_df_bed = pd.read_csv('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/phenotype_tracks/atac_coord_hg38.bed', sep='\t', header = None)
atac_bed=pybedtools.BedTool.from_dataframe(atac_df_bed)
coord_bed=pybedtools.BedTool.from_dataframe(coord_file)
intersect_bed=coord_bed.intersect(atac_bed,loj=True, wa = True, wb =True, s=False)
intersect_df=intersect_bed.to_dataframe(names=['chrom', 'start', 'end', 'name', 'score', 'strand','meta_id','chrom2', 'start2', 'end2', 'name_peak', 'score2', 'strand2',])
#%%
intersect_df=intersect_df[['chrom', 'start', 'end', 'name', 'score', 'strand','name_peak']]
#%%
intersect_df_info=intersect_df.merge(atac_df, left_on = ['name_peak'],right_on=['name'], how='inner')
# %%
metadata_sorted['metadata_sorted_order'] = metadata_sorted.index
intersect_metadata=intersect_df_info.merge(metadata_sorted, left_on='name_x', right_on='name', how='inner')
#%%
intersect_metadata_sorted=intersect_metadata.sort_values(['right_group','left_group','te_age','atac_group','ccres_group'])
inflame_curated_list = ['peak_33964','peak_33965','peak_88804','peak_96158','peak_495670','peak_301769','peak_23121','peak_391166','peak_354992','peak_142142','peak_304511','peak_304511','peak_228573','peak_275992','peak_601257','peak_168424','peak_486013','peak_168548','peak_60256','peak_93404','peak_371945','peak_187896','peak_133228','peak_200451','peak_334794','peak_329028','peak_183824']
intersect_metadata_sorted['inflamme_curated'] = 1
intersect_metadata_sorted.loc[intersect_metadata_sorted['name_peak'].isin(inflame_curated_list), 'inflamme_curated'] = 0
#%%
atac_enrich_matrix=intersect_metadata_sorted.iloc[:,13:106].to_numpy()
#%%
binding_indices=intersect_metadata_sorted[(intersect_metadata_sorted['right_group']==0)|(intersect_metadata_sorted['left_group']==0)].index
nonbinding_indices=intersect_metadata_sorted[(intersect_metadata_sorted['right_group']==1)&(intersect_metadata_sorted['left_group']==1)].index
only_right = intersect_metadata_sorted[(intersect_metadata_sorted['right_group']==0)&(intersect_metadata_sorted['left_group']==1)].index
only_left = intersect_metadata_sorted[(intersect_metadata_sorted['right_group']==1)&(intersect_metadata_sorted['left_group']==0)].index

intersect_right_left = intersect_metadata_sorted[(intersect_metadata_sorted['right_group']==0)&(intersect_metadata_sorted['left_group']==0)].index

#%%
intersect_metadata_sorted['motif_group'] = 0
intersect_metadata_sorted.loc[intersect_metadata_sorted.index.isin(intersect_right_left), 'motif_group'] = 0
intersect_metadata_sorted.loc[intersect_metadata_sorted.index.isin(only_right), 'motif_group'] = 1
intersect_metadata_sorted.loc[intersect_metadata_sorted.index.isin(only_left), 'motif_group'] = 2
intersect_metadata_sorted.loc[intersect_metadata_sorted.index.isin(nonbinding_indices), 'motif_group'] = 3
age_subgroups = np.unique(intersect_metadata_sorted['te_age'].sort_values())
age_subgroup = {subgroup: num for num, subgroup in enumerate(age_subgroups)}
age_anno=metadata_age['te_age'].map(age_subgroup)
subgroup_anno=intersect_metadata_sorted['motif_group']
subgroup_anno_atac = intersect_metadata_sorted['atac_group']
subgroup_anno_curated = intersect_metadata_sorted['inflamme_curated']
age_anno = intersect_metadata_sorted['te_age']
# %%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [atac_enrich_matrix,], 
    #alignment=alignment_sorted,
    #show_alignment=True, 
    heatmap_color=['viridis'],
    heatmap_mode='overlay', 
    vlim = [[0,100]], 
    annotation=True, 
    anno_col = ['Blues',['purple','red','blue','white'],['green','white'],['orange','white']], 
    annotation_data=[age_anno,subgroup_anno,subgroup_anno_curated,],
    anno_cbar_label=[age_subgroups,['both','right only','left only','no motif/others'],['curated','NA']],
    anno_cbar_title=['TEA-TIME','motif group','curated'], 
    opacity = 1.0,
    )
#%%
import scipy
distance_matrix = scipy.spatial.distance.pdist(atac_enrich_matrix, metric='euclidean')
linkage_array = scipy.cluster.hierarchy.linkage(distance_matrix, method = 'ward')
new_order = scipy.cluster.hierarchy.leaves_list(linkage_array)
atac_enrich_sorted=atac_enrich_matrix[new_order]
age_anno_sorted = age_anno[new_order]
subgroup_anno_sorted = subgroup_anno[new_order]
subgroup_anno_curated_sorted = subgroup_anno_curated[new_order]
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [atac_enrich_sorted,], 
    #alignment=alignment_sorted,
    #show_alignment=True, 
    heatmap_color=['viridis'],
    heatmap_mode='overlay', 
    vlim = [[0,100]], 
    annotation=True, 
    anno_col = ['Blues',['purple','red','blue','white'],['green','white'],['orange','white']], 
    annotation_data=[age_anno_sorted,subgroup_anno_sorted,subgroup_anno_curated_sorted,],
    anno_cbar_label=[age_subgroups,['both','right only','left only','no motif/others'],['curated','NA']],
    anno_cbar_title=['TEA-TIME','motif group','curated'], 
    opacity = 1.0,
    agg_major_tick=1,
    xlim=[50,85],
    xticklabels = intersect_metadata_sorted.iloc[:,13:106].columns,
    xticklabels_fs = 3,
    xticklabels_rt = 90
    )
#%%
from scipy.cluster.hierarchy import fcluster
cluster_numbers = fcluster(linkage_array, t=5, criterion='maxclust')
intersect_metadata_sorted['cluster_number'] = cluster_numbers
#%%
meta_clustered_sorted=intersect_metadata_sorted.sort_values(['motif_group','cluster_number'])
cluster_sort=meta_clustered_sorted.index
atac_enrich_sorted=atac_enrich_matrix[cluster_sort]
age_anno_sorted = age_anno[cluster_sort]
subgroup_anno_sorted = subgroup_anno[cluster_sort]
subgroup_anno_curated_sorted = subgroup_anno_curated[cluster_sort]
# %%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [atac_enrich_sorted,], 
    #alignment=alignment_sorted,
    #show_alignment=True, 
    heatmap_color=['viridis'],
    heatmap_mode='overlay', 
    vlim = [[0,100]], 
    annotation=True, 
    anno_col = ['Blues',['purple','red','blue','white'],['green','white'],['orange','white']], 
    annotation_data=[age_anno_sorted,subgroup_anno_sorted,subgroup_anno_curated_sorted,],
    anno_cbar_label=[age_subgroups,['both','right only','left only','no motif/others'],['curated','NA']],
    anno_cbar_title=['TEA-TIME','motif group','curated'], 
    opacity = 1.0,
    agg_major_tick=1,
    xlim=[50,85],
    xticklabels = intersect_metadata_sorted.iloc[:,13:106].columns,
    xticklabels_fs = 3,
    xticklabels_rt = 90
    )
#%%
import scipy
atac_zscore=scipy.stats.zscore(atac_enrich_matrix,axis=0)
distance_matrix = scipy.spatial.distance.pdist(atac_zscore, metric='euclidean')
linkage_array = scipy.cluster.hierarchy.linkage(distance_matrix, method = 'ward')
new_order = scipy.cluster.hierarchy.leaves_list(linkage_array)
atac_zscore_sorted=atac_zscore[new_order]
age_anno_sorted = age_anno[new_order]
subgroup_anno_sorted = subgroup_anno[new_order]
subgroup_anno_curated_sorted = subgroup_anno_curated[new_order]
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [atac_zscore_sorted,], 
    #alignment=alignment_sorted,
    #show_alignment=True, 
    heatmap_color=[custom_cmap.vlag_r_mpl],
    heatmap_mode='overlay', 
    vlim = [[-1,1]], 
    annotation=True, 
    anno_col = ['Blues',['purple','red','blue','white'],['green','white'],['orange','white']], 
    annotation_data=[age_anno_sorted,subgroup_anno_sorted,subgroup_anno_curated_sorted,],
    anno_cbar_label=[age_subgroups,['both','right only','left only','no motif/others'],['curated','NA']],
    anno_cbar_title=['TEA-TIME','motif group','curated'], 
    opacity = 1.0,
    agg_major_tick=1,
    colorbar=True,
    colorbar_steps=[0.1]
    #xlim=[50,85],
    #xticklabels = intersect_metadata_sorted.iloc[:,13:106].columns,
    #xticklabels_fs = 3,
    #xticklabels_rt = 90
    )
#%%
from scipy.cluster.hierarchy import fcluster
cluster_numbers = fcluster(linkage_array, t=25, criterion='maxclust')
intersect_metadata_sorted['cluster_number'] = cluster_numbers
#%%
meta_clustered_sorted=intersect_metadata_sorted.sort_values(['motif_group','cluster_number'])
cluster_sort=meta_clustered_sorted.index
atac_zscore_sorted=atac_zscore[cluster_sort]
age_anno_sorted = age_anno[cluster_sort]
subgroup_anno_sorted = subgroup_anno[cluster_sort]
subgroup_anno_curated_sorted = subgroup_anno_curated[cluster_sort]
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [atac_zscore_sorted,], 
    #alignment=alignment_sorted,
    #show_alignment=True, 
    heatmap_color=['bwr'],
    heatmap_mode='overlay', 
    vlim = [[-1,1]], 
    annotation=True, 
    anno_col = ['Blues',['purple','red','blue','white'],['green','white'],['orange','white']], 
    annotation_data=[age_anno_sorted,subgroup_anno_sorted,subgroup_anno_curated_sorted,],
    anno_cbar_label=[age_subgroups,['both','right only','left only','no motif/others'],['curated','NA']],
    anno_cbar_title=['TEA-TIME','motif group','curated'], 
    opacity = 1.0,
    agg_major_tick=10,
    #xlim=[50,85],
    #xticklabels = intersect_metadata_sorted.iloc[:,13:106].columns,
    #xticklabels_fs = 3,
    #xticklabels_rt = 90
    )
# %%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [atac_enrich_matrix,], 
    #alignment=alignment_sorted,
    #show_alignment=True, 
    heatmap_color=['viridis'],
    heatmap_mode='overlay', 
    vlim = [[0,100]], 
    annotation=True, 
    anno_col = ['Blues',['purple','red','blue','white'],['green','white'],['orange','white']], 
    annotation_data=[age_anno,subgroup_anno,subgroup_anno_curated,subgroup_anno_atac],
    anno_cbar_label=[age_subgroups,['both','right only','left only','no motif/others'],['inflamme action','NA'],['ATAC peak','no ATAC']],
    anno_cbar_title=['TEA-TIME','motif group','curated','ATAC peak'], 
    opacity = 1.0,

    )
