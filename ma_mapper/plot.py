import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap
from matplotlib.pyplot import gcf
from matplotlib.colors import BoundaryNorm
import os
import sys
import numpy as np
import pandas as pd
from . import mapper
from . import logger
if sys.version_info >= (3, 8, 0):
    from typing import Literal
else:
    from typing_extensions import Literal
#%%
_ACOLARG = Literal['dna', 'white']
def overlay_plot(list_of_overlays: list, 
                 save_to_file: bool = False, 
                 nucleotide_color:_ACOLARG = 'dna',
                 output_filepath: str | None = None, 
                 image_res: int = 600,
                 metadata: str | pd.DataFrame | None = None,
                 age_annotation_cmap: str = 'Blues',
                 plot_title: str = 'alignment_plot',
                 **kwargs) -> None:
    plt.rcParams['figure.dpi'] = image_res
    
    fig_width = 8
    fig_height = 5
    grid_nrow = 10
    grid_ncol = 25 
    #grid parameters = row_start, row_end, col_start, col_end
    mainplot_grid_params = (0,10,0,22)
    align_legend_params = (0,3,23,24)

    fig = plt.figure(figsize=(fig_width,fig_height))
    grid = fig.add_gridspec(nrows = grid_nrow, ncols = grid_ncol, hspace=0)
    heatmap_grid = fig.add_subplot(grid[mainplot_grid_params[0]:mainplot_grid_params[1],mainplot_grid_params[2]:mainplot_grid_params[3]])
    heatmap_grid.set_yticks([])
    heatmap_grid.set_title(plot_title)

    if nucleotide_color == 'dna':
        nucleotide_labels = ['gap', 'A', 'C', 'T', 'G']
        nucleotide_color_list = ['grey','green','yellow','red','blue']
        align_cmap = ListedColormap(nucleotide_color_list)
    elif nucleotide_color == 'white':
        align_cmap =ListedColormap(['grey','white','white','white','white'])

    for idx, overlay in enumerate(list_of_overlays):
        if idx == 0:
            heatmap_core = heatmap_grid.imshow(overlay, aspect = 'auto',cmap= align_cmap, interpolation='nearest', vmax=5)
        else:
            heatmap_grid.imshow(overlay, aspect = 'auto')

    if nucleotide_color == 'dna':
        align_legend_grid = fig.add_subplot(grid[align_legend_params[0]:align_legend_params[1],align_legend_params[2]:align_legend_params[3]])
        align_legend=fig.colorbar(heatmap_core, cax=align_legend_grid, ticks=[0.5, 1.5, 2.5, 3.5, 4.5], orientation='vertical')
        align_legend.set_ticklabels(nucleotide_labels, fontsize='small')
        align_legend.ax.set_title('nucleotide', fontsize='small')

    if metadata is not None:
        if metadata is str:
            metadata=pd.read_csv(metadata)
        if 'te_age' not in metadata:
            logger.warning('cannot find te_age in metadata, trying to find from reference_table')
            reference_tbl = kwargs.get('reference_table', None)
            metadata_age = mapper.match_age_to_id_metadata(metadata, reference_tbl)
        else:
            metadata_age = metadata

        age_annot_grid_params = (0,10,21,22)
        age_annot_legend_params = (5,10,23,24)

        age_unique=metadata_age.te_age.sort_values().unique()
        age_boundary = np.concatenate(([-0.5], age_unique))
        age_cmap= plt.get_cmap('Blues', len(age_unique))
        norm = BoundaryNorm(age_boundary, len(age_unique))

        age_annot_grid = fig.add_subplot(grid[age_annot_grid_params[0]:age_annot_grid_params[1],age_annot_grid_params[2]:age_annot_grid_params[3]])
        age_annot = age_annot_grid.imshow(metadata_age.te_age.values.reshape(-1, 1), aspect = 'auto',cmap= age_cmap, interpolation='nearest', norm=norm)
        age_annot_grid.set_xticks([])
        age_annot_grid.set_yticks([])
        
        age_annot_legend_grid = fig.add_subplot(grid[age_annot_legend_params[0]:age_annot_legend_params[1],age_annot_legend_params[2]:age_annot_legend_params[3]])
        age_annot_legend = fig.colorbar(age_annot, cax = age_annot_legend_grid)
        age_annot_legend.set_ticks([(age_boundary[i] + age_boundary[i+1]) / 2 for i in range(len(age_boundary) - 1)])
        age_annot_legend.set_ticklabels(age_unique, fontsize = 'small')
        age_annot_legend.ax.set_title('te_age', fontsize ='small')

    if save_to_file == False:
        plt.rcParams['savefig.dpi'] = image_res
        plt.show()
    else: 
        if output_filepath is None:
            output_filepath = os.path.dirname(os.path.abspath(__file__)) + '/overlay.png'
        fig.write_image(output_filepath)
