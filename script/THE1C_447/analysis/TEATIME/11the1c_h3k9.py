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
# %%
bam_file =  '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/H3K9-bam/H3K9me3_HEK293_rep1.bam'
#%%
bam_forward=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_forward', pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10, offset = 100)
bam_forward=bam_forward[noNA_indices]
bam_reverse=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_reverse', pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10, offset = 100)
bam_reverse=bam_reverse[noNA_indices]
mean_forward=mapper.normalise(alignment=alignment_filtered, mapped_data=bam_forward)
mean_reverse=mapper.normalise(alignment=alignment_filtered, mapped_data=bam_reverse)
#%%
bam_min=np.minimum(mean_forward,mean_reverse)
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
    vlim = [[0,0.00000001]], 
    hm_opacity = 0.9, 
    hm_transparency_mode = 'gradient', 
    aggregated=True, 
    aggregated_data=[bam_min], 
    agg_colset=['purple'],
    ag_opacity=0.5,
    agg_ylim=[[0,0.01]],
    agg_titles=['H3K3 ChIP-seq','znf808 ChIP-exo reverse'], 
    agg_ylabel=['signal coverage','signal coverage'],
    agg_xlabel = 'position (bp)',
    #annotation=True, 
    hm_plot_title = 'H3K9 ChIP-seq signals\nfrom HEK239 cell lines on THE1C MSA',
    hm_xlabel = 'position (bp)',
    hm_ylabel = 'seqeunces',
    #anno_col = ['Blues'], 
    #annotation_data=[age_anno],
    #anno_cbar_label=[age_subgroups],
    #anno_cbar_title=['TEA-TIME'], 
    #colorbar=False,
    agg_major_tick=50
    )
#%%