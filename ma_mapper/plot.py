import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm
import matplotlib
import os
import sys
import numpy as np
import pandas as pd
from . import mapper
from . import logger
if sys.version_info >= (3, 8, 0):
    from typing import Literal, Tuple, List
else:
    from typing_extensions import Literal, Tuple, List
#%%
_ACOLARG = Literal['dna', 'white']
def overlay_heatmap_handler(heatmap_grid, overlay, cmap, idx, data_vmax,data_vmin):
    
    if data_vmax is not None and data_vmin is not None:
        heatmap_grid.imshow(overlay, aspect='auto', cmap=cmap, vmax=data_vmax[idx-1], vmin=data_vmin[idx-1])
    elif data_vmax is not None:
        heatmap_grid.imshow(overlay, aspect='auto', cmap=cmap, vmax=data_vmax[idx-1])
    elif data_vmin is not None:
        heatmap_grid.imshow(overlay, aspect='auto', cmap=cmap, vmin=data_vmin[idx-1])
    else:
        heatmap_grid.imshow(overlay, aspect='auto', cmap=cmap)


def overlay_plot(list_of_overlays: list, 
                 save_to_file: bool = False, 
                 nucleotide_color:_ACOLARG = 'dna',
                 output_filepath: str | None = None, 
                 image_res: int = 600,
                 metadata: str | pd.DataFrame | None = None,
                 age_annotation:bool = True,
                 age_annotation_cmap: str = 'Blues',
                 data_cmap: List|None = None,
                 data_vmin: List|None = None,
                 data_vmax: List|None = None,
                 plot_title: str = 'alignment_plot',
                 show_nucleotide_legend:bool = False,
                 show_data_legend:bool = False,
                 data_label: List|None = None,
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
            heatmap_core = heatmap_grid.imshow(overlay, aspect = 'auto',cmap= align_cmap, interpolation='nearest', vmin=0,vmax=5)
        else:
            ncolors = 256
            if data_cmap is not None:
                overlay_color = data_cmap[idx-1]
            else:
                overlay_color = 'Blues'
            overlay_cmap = plt.get_cmap(overlay_color)(range(ncolors))
            overlay_cmap[:,-1] = np.linspace(0,1.0,ncolors)
            overlay_cmap_alpha = ListedColormap(colors=overlay_cmap)
            overlay_heatmap_handler(heatmap_grid, overlay, overlay_cmap_alpha, idx, data_vmax, data_vmin)
            
    if show_data_legend == True:
        show_nucleotide_legend = False #data legend should always override 
        number_of_overlay = (len(list_of_overlays)-1)
        if number_of_overlay == 1: #set minimum overlay number for gridspec
            number_of_overlay = 2
        align_legend_grid = fig.add_subplot(grid[align_legend_params[0]:int(number_of_overlay/2),align_legend_params[2]:align_legend_params[3]])
        data_pallete = []
        for cmap_name in data_cmap:
            cmap=plt.get_cmap(cmap_name)
            data_pallete.append(cmap(cmap.N))
        overlay_cmap = ListedColormap(colors=data_pallete)
        bounds = range(overlay_cmap.N+1)
        norm = BoundaryNorm(bounds, overlay_cmap.N)
        data_legend=fig.colorbar(matplotlib.cm.ScalarMappable(cmap=overlay_cmap, norm=norm), cax=align_legend_grid, orientation='vertical')
        data_legend.set_ticks([(bounds[i] + bounds[i+1]) / 2 for i in range(len(bounds) - 1)])
        if data_label is None:
            data_label = []
            for i in range(len(list_of_overlays)-1):
                data_label.append('layer_'+ str(i))
        data_legend.set_ticklabels(data_label, fontsize = 'small')
        data_legend.ax.set_title('data', fontsize ='small')

    if show_nucleotide_legend == True:
        if nucleotide_color == 'dna':
            align_legend_grid = fig.add_subplot(grid[align_legend_params[0]:align_legend_params[1],align_legend_params[2]:align_legend_params[3]])
            align_legend=fig.colorbar(heatmap_core, cax=align_legend_grid, ticks=[0.5, 1.5, 2.5, 3.5, 4.5], orientation='vertical')
            align_legend.set_ticklabels(nucleotide_labels, fontsize='small')
            align_legend.ax.set_title('nucleotide', fontsize='small')



    if metadata is not None:
        if metadata is str:
            metadata=pd.read_csv(metadata)
        if age_annotation:
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
                age_cmap= plt.get_cmap(age_annotation_cmap, len(age_unique))
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
