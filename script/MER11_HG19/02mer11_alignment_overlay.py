#%%
import sys
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
subfamily = 'MER11'
from ma_mapper import mapper
alignment_file = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11_HG19/mer11_w_peak.fasta.aligned'
import pandas as pd
from ma_mapper import sequence_alignment
#%% load TE coords
hg19_coord_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11_HG19/t28_final_peaks_te_mer11_coords.bed'
hg19_coord_df=pd.read_csv(hg19_coord_file, sep='\t', header=None)
hg19_coord_df.columns=['chrom','start','end','subfamily','#peak']
hg19_coord_df['start'] = hg19_coord_df['start'] -1
#%% load repeatmasker table
repeatmasker_hg19 = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/repeatmasker/repeatmasker/hg19_repeatlib2014/hg19.fa.out.tsv'
repeatmasker_hg19_df=pd.read_csv(repeatmasker_hg19, sep='\t', index_col=0)
#%% make bed file for seqeuence extraction and alignment
prealign_metadata=hg19_coord_df.merge(repeatmasker_hg19_df, how='left',left_on=['chrom','start','end'], right_on=['genoName','genoStart','genoEnd'])[['chrom','start','end','subfamily','#peak','strand']]
prealign_metadata['name'] = 'MER11_'+prealign_metadata.index.astype(str)
prealign_metadata['score'] = 10
coord_file = prealign_metadata[['chrom','start','end','name','score','strand']]
#%% parse and filter alignment into matrix
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file, col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
metadata_filtered
#%% merge prealignmetadata to add more information -> make plot annotation
metadata_w_info=metadata_filtered.merge(prealign_metadata[['subfamily','name','#peak']], how='left', on='name')

#%% make annotation
import numpy as np
subgroups = np.unique(metadata_w_info['subfamily'].astype(str))
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno=metadata_w_info['subfamily'].map(numerical_subgroup)
#%% plot
from ma_mapper import plots
import importlib
importlib.reload(plots)
plots.plot_experimental(
    show_alignment=True,
    alignment=alignment_filtered, 
    alignment_col='dna', 
    show_alignment_colbar=True,
    colorbar=True,  
    agg_major_tick=100, 
    annotation=True, 
    annotation_data=[subgroup_anno], 
    anno_col=[['red','yellow','blue']], 
    anno_title=['subfamily'],
    anno_cbar=True, 
    anno_cbar_label=[['MER11A','MER11B','MER11C']], )
#%%
bed_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/kzfp_peak_bed/hg19_kzfps_combined.bed'
kzfp_df=pd.read_csv(bed_filepath, sep='\t', header=None)
kzfp_df.columns=['chrom','start','end','name','score','strand']
#%%
#%%

#%%
bed_file=kzfp_df[kzfp_df['name'].str.contains('ZNF808')]
znf808=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=False, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
bed_file=kzfp_df[kzfp_df['name'].str.contains('ZNF525')]
znf525=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=False, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
bed_file=kzfp_df[kzfp_df['name'].str.contains('ZNF727')]
znf727=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=False, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
bed_file=kzfp_df[kzfp_df['name'].str.contains('ZNF468')]
znf468=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=False, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
bed_file=kzfp_df[kzfp_df['name'].str.contains('ZNF440')]
znf440=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=False, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
bed_file=kzfp_df[kzfp_df['name'].str.contains('ZNF433')]
znf433=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=False, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
bed_file=kzfp_df[kzfp_df['name'].str.contains('ZNF611')]
znf611=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=False, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
bed_file=kzfp_df[kzfp_df['name'].str.contains('ZNF578')]
znf578=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=False, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
#%%
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
def create_transparent_colormap(base_color, alpha=0.8):
    # Define the RGBA colors - white with full transparency and the base color with specified transparency
    colors = [(1, 1, 1, 0), plt.cm.colors.to_rgba(base_color, alpha)]
    # Create the ListedColormap
    cmap = ListedColormap(colors)
    return cmap
znf808_col=create_transparent_colormap('blue')
znf525_col=create_transparent_colormap('green')
znf727_col=create_transparent_colormap('darkviolet')
znf468_col=create_transparent_colormap('dodgerblue')
znf440_col=create_transparent_colormap('gold')
znf433_col=create_transparent_colormap('yellowgreen')
znf611_col=create_transparent_colormap('orangered')
znf578_col=create_transparent_colormap('red')

# %%
importlib.reload(plots)
plots.plot_experimental(
    data = [znf808,znf525,znf727,znf468,znf440,znf433,znf611,znf578], 
    heatmap_color=[znf808_col,znf525_col,znf727_col,znf468_col,znf440_col,znf433_col,znf611_col,znf578_col,],
    heatmap_mode='overlay', 
    vlim = [[0,7],[0,7],[0,7],[0,7],[0,7],[0,7],[0,7],[0,7]], 
    opacity = 1.0, 

    show_alignment=True,
    alignment=alignment_filtered, 
    alignment_col='nulc_white', 

    agg_major_tick=100, 
    annotation=True, 
    annotation_data=[subgroup_anno], 
    anno_col=[['red','yellow','blue']], 
    anno_title=['subfamily'],
    anno_cbar=True, 
    anno_cbar_label=[['MER11A','MER11B','MER11C']], )
# %%
target_order=pd.read_csv('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11_HG19/t28_final_peaks_te_mer11_counts_kzfps.csv')
# %%
merged_df=target_order.merge(metadata_w_info, how='left', left_on='name', right_on='#peak')
# %%
sorted_index=merged_df['original_order']
# %%
znf808_sorted = znf808[sorted_index]
znf525_sorted = znf525[sorted_index]
znf727_sorted = znf727[sorted_index]
znf468_sorted = znf468[sorted_index]
znf440_sorted = znf440[sorted_index]
znf433_sorted = znf433[sorted_index]
znf611_sorted = znf611[sorted_index]
znf578_sorted = znf578[sorted_index]
alignment_sorted = alignment_filtered[sorted_index]
metadata_sorted = metadata_w_info.reindex(sorted_index)
subgroup_anno_sorted = metadata_sorted['subfamily'].map(numerical_subgroup)
# %%
importlib.reload(plots)
plots.plot_experimental(
    data = [znf808_sorted,znf525_sorted,znf727_sorted,znf468_sorted,znf440_sorted,znf433_sorted,znf611_sorted,znf578_sorted], 
    heatmap_color=[znf808_col,znf525_col,znf727_col,znf468_col,znf440_col,znf433_col,znf611_col,znf578_col,],
    heatmap_mode='overlay', 
    vlim = [[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1]], 
    opacity = 1.0, 
    
    show_alignment=True,
    alignment=alignment_sorted, 
    alignment_col='nulc_white', 
    ylim = [100,0],
    agg_major_tick=100, 
    annotation=True, 
    annotation_data=[subgroup_anno_sorted], 
    anno_col=[['red','yellow','blue']], 
    anno_title=['subfamily'],
    anno_cbar=True, 
    anno_cbar_label=[['MER11A','MER11B','MER11C']], 
    )
# %%
np.savetxt("/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11_HG19/delivery_01/alignment_sorted.txt", alignment_sorted, delimiter="\t")
np.savetxt("/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11_HG19/delivery_01/znf808_sorted.txt", znf808_sorted, delimiter="\t")
np.savetxt("/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11_HG19/delivery_01/znf525_sorted.txt", znf525_sorted, delimiter="\t")
np.savetxt("/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11_HG19/delivery_01/znf727_sorted.txt", znf727_sorted, delimiter="\t")
np.savetxt("/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11_HG19/delivery_01/znf468_sorted.txt", znf468_sorted, delimiter="\t")
np.savetxt("/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11_HG19/delivery_01/znf440_sorted.txt", znf440_sorted, delimiter="\t")
np.savetxt("/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11_HG19/delivery_01/znf433_sorted.txt", znf433_sorted, delimiter="\t")
np.savetxt("/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11_HG19/delivery_01/znf611_sorted.txt", znf611_sorted, delimiter="\t")
np.savetxt("/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11_HG19/delivery_01/znf578_sorted.txt", znf578_sorted, delimiter="\t")
#%%
metadata_sorted.to_csv('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11_HG19/delivery_01/metadata_sorted.txt',index=False,sep='\t')
#%%