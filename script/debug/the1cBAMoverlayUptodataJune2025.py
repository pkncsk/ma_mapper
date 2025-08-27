#%%
from networkx import radius
from ma_mapper import mapper
from ma_mapper import plots
from ma_mapper import custom_cmap
import numpy as np
import pandas as pd
#%%
subfamily='THE1C'
alignment_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/alignment/THE1C.fasta.aligned'
genomewide_data_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/sequenceing_data_fastqsambam/znf267.sorted.bam'
#%%
filtered_alignment_matrix, alignment_coordinate  = mapper.parse_and_filter(alignment_file=alignment_filepath,col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
data_matrix_forward=mapper.map_and_overlay(alignment_filepath, genomewide_data_filepath,data_format='read_forward', col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
data_matrix_reverse=mapper.map_and_overlay(alignment_filepath, genomewide_data_filepath,data_format='read_reverse', col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
#%%
from scipy.ndimage import gaussian_filter1d
# Apply Gaussian blur row-wise
sigma = 5  # standard deviation of the Gaussian kernel
rad=5
blurred_forward = gaussian_filter1d(data_matrix_forward, sigma=sigma, axis=0, radius=rad)
blurred_reverse = gaussian_filter1d(data_matrix_reverse, sigma=sigma, axis=0, radius=rad)
#%%
from matplotlib import cm
from matplotlib.colors import ListedColormap
brighter_cmap = cm.get_cmap('RdBu')(np.linspace(0.1, 0.9, 256))
from ma_mapper import plots
from ma_mapper import mapper
plots.plot(
    data = [blurred_forward,blurred_reverse], 
    alignment=filtered_alignment_matrix,
    show_alignment=False, 
    heatmap_color=['Reds','Blues'],
    heatmap_mode='overlay',
    heatmap_title=['ZNF267 ChIP-exo on THE1C MSA'],
    heatmap_title_fs = 10, 
    heatmap_xlabel = 'position (bp)',
    heatmap_xlabel_fs = 10,
    heatmap_ylabel = 'sequences',
    heatmap_ylabel_fs = 10,
    vlim = [[0,0.125],[0,0.125]], 
    opacity = 0.9, 
    hm_transparency_mode= 'constant',
    hm_interpolation = 'auto',
    colorbar=True,
    colorbar_steps=[0.025,0.025],
    save_to_file='bam_blurred.png',
    image_res=300
    )
#%%
from ma_mapper import plots
from ma_mapper import mapper
plots.plot(
    data = [data_matrix_forward,data_matrix_reverse], 
    alignment=filtered_alignment_matrix,
    show_alignment=False, 
    heatmap_color=['Reds','Blues'],
    heatmap_mode='overlay',
    heatmap_title=['ZNF267 ChIP-exo on THE1C MSA'],
    heatmap_title_fs = 10, 
    heatmap_xlabel = 'position (bp)',
    heatmap_xlabel_fs = 10,
    heatmap_ylabel = 'sequences',
    heatmap_ylabel_fs = 10,
    vlim = [[0,0.125],[0,0.125]], 
    opacity = 0.9, 
    hm_transparency_mode= 'constant',
    #hm_interpolation = 'auto',
    colorbar=True,
    colorbar_steps=[0.025,0.025]
    
    )
# %%
