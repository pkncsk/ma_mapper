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
    data = [ap1,atac], 
    alignment=alignment_filtered,
    show_alignment=True, 
    heatmap_color=['Greens','Purples'],
    heatmap_mode='overlay', 
    vlim = [[0,7]], 
    opacity = 0.5, 
    aggregated=True, 
    aggregated_data=[cov_atac,cov_ap1], 
    agg_colset=['purple','green'],
    agg_ylim=[[0,70],[0,70]],
    agg_titles=['AP-1 motif','znf267-p65 motif'], 
    agg_ylabel=['perc_coverage', 'perc_coverage'],
    annotation=True, 
    anno_col = ['Blues'], 
    annotation_data=[age_anno],
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
metadata_age['right_group'] = 0
metadata_age.loc[metadata_age.index.isin(binding_indices_right), 'right_group'] = 1
metadata_age.loc[metadata_age.index.isin(nonbinding_indices_right), 'right_group'] = 0
metadata_age['left_group'] = 0
metadata_age.loc[metadata_age.index.isin(binding_indices_left), 'left_group'] = 1
metadata_age.loc[metadata_age.index.isin(nonbinding_indices_left), 'left_group'] = 0
#%%
metadata_age['motif_group'] = 0
metadata_age.loc[(metadata_age['right_group'] == 0) & (metadata_age['left_group'] == 0), 'motif_group'] = 0
metadata_age.loc[(metadata_age['right_group'] == 0) & (metadata_age['left_group'] == 1), 'motif_group'] = 1
metadata_age.loc[(metadata_age['right_group'] == 1) & (metadata_age['left_group'] == 0), 'motif_group'] = 2
metadata_age.loc[(metadata_age['right_group'] == 1) & (metadata_age['left_group'] == 1), 'motif_group'] = 3

#%%
import scipy
peaks, _ = scipy.signal.find_peaks(cov_atac, width = 6)
# %%
annot_indices = np.unique(np.where(atac[:, peaks] != 0)[0])
nonannot_indices=list(set(np.arange(atac.shape[0])) - set(annot_indices))
#%%
metadata_age['atac_group'] = 0
metadata_age.loc[metadata_age.index.isin(annot_indices), 'atac_group'] = 1
metadata_age.loc[metadata_age.index.isin(nonannot_indices), 'atac_group'] = 0
atac_sort=metadata_age.sort_values(['atac_group']).index
# %%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [atac[atac_sort]], 
    alignment=alignment_filtered,
    show_alignment=True, 
    heatmap_color=['viridis'],
    heatmap_mode='overlay', 
    vlim = [[0,1]], 
    opacity = 1.0, 
    annotation=True, 
    anno_col = ['Blues'], 
    annotation_data=[age_anno[atac_sort]],
    anno_cbar_label=[age_subgroups],
    anno_cbar_title=['TEA-TIME'], 
    colorbar=False,
    )
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
atac_the1_df_lean = atac_the1_df[['name','TSS','TSS_5000']].rename(columns={'name':'name_peak'})
# %%
intersect_info=intersect_df.merge(atac_the1_df_lean, on='name_peak')
# %%
metadata_info=metadata_age.merge(intersect_info, on=['chrom','start','end','name','score','strand'], how='left')
# %%
metadata_sorted=metadata_info.sort_values(['motif_group','te_age','atac_group'], ascending=False)
sorted_indices=metadata_sorted.index
age_subgroups = np.unique(metadata_sorted['te_age'].sort_values())
age_subgroup = {subgroup: num for num, subgroup in enumerate(age_subgroups)}
age_anno=metadata_sorted['te_age'].map(age_subgroup)
motif_anno = metadata_sorted['motif_group']
atac_anno =metadata_sorted['atac_group']
tss_anno = metadata_sorted['TSS']
tss5k_anno = metadata_sorted['TSS_5000']
phylop_sorted=phylop_447[sorted_indices]
alignment_sorted=alignment_filtered[sorted_indices]
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [phylop_sorted,], 
    alignment=alignment_sorted,
    #show_alignment=True, 
    heatmap_color=['RdBu'],
    heatmap_mode='overlay', 
    vlim = [[-1.5,1.5]], 
    opacity = 1.0, 
    annotation=True, 
    anno_col = ['Blues',['white','yellow','green','orange'],['white','blue'],['white','red'],['white','purple']], 
    annotation_data=[age_anno,motif_anno,atac_anno,tss_anno,tss5k_anno],
    anno_cbar_label=[age_subgroups,['no motif','left only','right only','both'],['NA','ATAC'],['NA','TSS'],['NA','TSS_5k']],
    anno_cbar_title=['TEA-TIME','motif group', 'ATAC','TSS','TSS_5k'], 
    colorbar=True,
    colorbar_steps=[0.1],
    agg_xlabel='consensus position (bp)'
    )
# %%
metadata_age_filtered=metadata_sorted
binding_indices = metadata_age_filtered[(metadata_age_filtered['motif_group'].isin([1,2,3]))&(metadata_age_filtered['atac_group']==1)].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['motif_group']==0)&(metadata_age_filtered['atac_group']==0)].index
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
    agg_ylim=[[-1.5,1.5],[0,15],[-15,0]],
    agg_yhighlight=[[98,108],[304,314]],
    agg_yhighlight_col= ['green','green'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    agg_xlabel='consensus position (bp)'
    #xlim=[304,313],
    #gg_major_tick=1,
    )
# %%
metadata_age_filtered=metadata_sorted
binding_indices = metadata_age_filtered[(metadata_age_filtered['motif_group'].isin([1,2,3]))&(metadata_age_filtered['atac_group']==1)].index
nonbinding_indices = metadata_age_filtered[(metadata_age_filtered['motif_group']==0)&(metadata_age_filtered['atac_group']==0)].index
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
    agg_colset=['grey','blue','red'],
    agg_ylim=[[-1.5,1.5],[0,15],[-15,0]],
    agg_yhighlight=[[98,108],[304,314]],
    agg_yhighlight_col= ['green','green'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_titles=['mean_phyloP','conserved\nmotif>no motif','accelerated\nmotif>no motif'], 
    agg_ylabel=[None,'-log10p_adj','log10p_adj'],
    colorbar=False,
    agg_xlabel='consensus position (bp)'
    #xlim=[304,313],
    #gg_major_tick=1,
    )
#%%