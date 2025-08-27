from configparser import Interpolation
from ma_mapper import mapper
from ma_mapper import plots
from ma_mapper import custom_cmap
import numpy as np
import pandas as pd
#%%
subfamily='THE1C'
alignment_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/alignment/THE1C.fasta.aligned'
genomewide_data_dir  = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447'
#%%
filtered_alignment_matrix, alignment_coordinate  = mapper.parse_and_filter(alignment_file=alignment_filepath,col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10,)
data_matrix=mapper.map_and_overlay(alignment_filepath, genomewide_data_dir,data_format='maf', count_arg='ref_freq', separated_maf=True,target_species='hg38', col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
#%%
from matplotlib import cm
from matplotlib.colors import ListedColormap

brighter_cmap = cm.get_cmap('RdBu')(np.linspace(0.05, 0.95, 256))
from ma_mapper import plots
from ma_mapper import custom_cmap
plots.plot(
    data = [data_matrix], 
    alignment=filtered_alignment_matrix,
    show_alignment=False, 
    heatmap_color=['viridis'],
    heatmap_mode='overlay',
    heatmap_title=['Common base frequency\nfrom Zoonomia dataset on THE1C MSA'],
    heatmap_title_fs = 10, 
    heatmap_xlabel = 'position (bp)',
    heatmap_xlabel_fs = 10,
    heatmap_ylabel = 'sequences',
    heatmap_ylabel_fs = 10,
    vlim = [[0,1.0]], 
    opacity = 0.999, 
    hm_transparency_mode= 'static',
    colorbar=True,
    colorbar_steps=[0.1],
    image_res=300,
    save_to_file='./MAF'
    )
#%%