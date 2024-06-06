#%%
import sys
from tabnanny import verbose
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import mapper
import pandas as pd
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
if sys.version_info >= (3, 8, 0):
    from typing import Literal, Tuple, List
else:
    from typing_extensions import Literal, Tuple, List
#%%


#%%
def plot_heatmap(data: list|np.ndarray|None = None,
                 save_to_file: bool = False,
                 output_filepath: str | None = None, 
                 image_res: int = 600, 
                 cmap: ListedColormap|str|None = None,
                 plot_title: str|None = None,
                 xlim: list|None = None,
                 ylim: list|None = None,
                 figsize:list=[2,2],
                 opacity: float = 1.0,
                 show_plot:bool = False, 
                 matplot_axes = None,
                 zero_col_centering:bool=False,
                 **kwargs):
    if data is None or not isinstance(data, (list, np.ndarray)):
        raise ValueError("Invalid data")
    
    cmap_dict = {'nulc_white': ['grey', 'white'], 'dna': ['grey', 'green', 'yellow', 'red', 'blue']}

    if isinstance(cmap, ListedColormap|LinearSegmentedColormap):
        plot_cmap_a = cmap
    elif cmap in cmap_dict:
        zero_col_centering=False
        plot_cmap = ListedColormap(cmap_dict[cmap])
        kwargs['vmax'] = 5
        kwargs['vmin'] = 0
        plot_cmap_a = plot_cmap
    else:
        cmap = cmap or 'Blues' #assign Blues to nonetype 
        ncolors = 256
        plot_cmap = plt.get_cmap(cmap)(range(ncolors))
        plot_cmap[:,-1] = np.linspace(0, opacity, ncolors)
        plot_cmap_a = ListedColormap(colors=plot_cmap)
    if zero_col_centering:
        vmax = kwargs.get('vmax', np.nanmax(data))
        vmin = kwargs.get('vmin', np.nanmin(data))
        norm = matplotlib.colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
        kwargs['norm'] = norm
        kwargs.pop('vmax', None)
        kwargs.pop('vmin', None)

    ax = matplot_axes if matplot_axes is not None else plt.subplots(figsize=figsize)[1]
    heatmap = ax.imshow(data, aspect='auto', cmap=plot_cmap_a, **kwargs)
    if plot_title is not None:
        ax.set_title(plot_title)
    if xlim is not None:
        ax.set_xlim(xlim[0],xlim[1])
    if ylim is not None:
        ax.set_ylim(ylim[0],ylim[1])
    if show_plot == True:
        plt.rcParams['figure.dpi'] = image_res
        plt.show()
    if save_to_file == True:
        plt.rcParams['savefig.dpi'] = image_res
        if output_filepath is None:
            output_filepath = f'{os.path.dirname(os.path.abspath(__file__))}/overlay.png'
        plt.savefig(output_filepath) 
# %%
_DATAMODE = Literal['1D','average']
def plot_bar(data: list|np.ndarray|None = None,
             alignment: np.ndarray|None = None,
             save_to_file: bool = False,
             output_filepath: str | None = None, 
             image_res: int = 600, 
             color: str|None = None,
             bar_title: str|None = None,
             xlim: list|None = None,
             ylim: list|None = None,
             figsize:list=[2,2],
             opacity: float = 1.0,
             show_plot:bool = False, 
             mode:_DATAMODE = 'average',
             matplot_axes = None,
             **kwargs):
    if data is None or not isinstance(data, (list, np.ndarray)):
        raise ValueError("Invalid data")
    if mode == '1D':
        dataplot = data
    elif mode == 'average':
        if alignment is None:
            raise ValueError("need alignment when mode is 'average'")
        nm_kwargs = {key: value for key, value in kwargs.items() if key in ('method')}
        dataplot=mapper.normalise(alignment=alignment,mapped_data=data,**nm_kwargs)
    if color is None:
        color = 'grey'
    ax = matplot_axes if matplot_axes is not None else plt.subplots(figsize=figsize)[1]
    bar = ax.bar(range(len(dataplot)), dataplot, color = color,alpha =opacity, **kwargs)
    ax.margins(x=0, y=0)
    if bar_title is not None:
        ax.set_title(bar_title)
    if xlim is not None:
        ax.set_xlim(xlim[0],xlim[1])
    if ylim is not None:
        ax.set_ylim(ylim[0],ylim[1])
    if show_plot == True:
        plt.rcParams['figure.dpi'] = image_res
        plt.show()
    if save_to_file == True:
        plt.rcParams['savefig.dpi'] = image_res
        if output_filepath is None:
            output_filepath = f'{os.path.dirname(os.path.abspath(__file__))}/barplot.png'
        plt.savefig(output_filepath)  
#%%
_ORIENT = Literal['horizontal','vertical']
def plot_colorbar(cmap: str|ListedColormap|LinearSegmentedColormap,
                  vmin:int|float|None = None,
                  vmax:int|float|None = None,
                  data:list|np.ndarray|None = None,
                  save_to_file: bool = False,
                  output_filepath: str | None = None, 
                  image_res: int = 600, 
                  cbar_label: str|None = None,
                  cbar_title: str|None = None,
                  tick_label: list|None = None,
                  figsize:list=[2,2],
                  show_plot:bool = False, 
                  centering:bool = False,
                  step: int|float|None = None,
                  orientation:_ORIENT = 'vertical',
                  matplot_axes = None,
                  **kwargs):
    import matplotlib as mpl
    cmap_dict = {'nulc_white': ['grey', 'white'], 'dna': ['grey', 'green', 'yellow', 'red', 'blue'],'dna_jaspar': ['white', 'green', 'blue', 'red', 'yellow']}
    if cmap is None:
        raise ValueError('need colormap')
    elif isinstance(cmap, ListedColormap|LinearSegmentedColormap):
        colorbar_cmap = cmap
    
    elif cmap in cmap_dict:
        colorbar_cmap = ListedColormap(cmap_dict[cmap])
        vmax = int(5)
        vmin = int(0)
        kwargs['ticks'] = [0.5, 1.5, 2.5, 3.5, 4.5]
    else:
        ncolors = 256
        plot_cmap = plt.get_cmap(cmap)(range(ncolors))
        colorbar_cmap = ListedColormap(colors=plot_cmap)

    if vmin is None and vmax is None:
        if data is not None:
            vmax = np.nanmax(data)
            vmin = np.nanmin(data)
        else:
            raise ValueError('insufficient input, please provide either vmin+vmax or data')
    if centering:
        norm = mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    else:
        if step is None:
            if (vmax is int and vmin is int) or cmap in cmap_dict:
                step =1
            else:
                step =0.1
        print(step)
        bounds =np.arange(vmin,vmax+step,step)
        extend_arg = 'neither'
        if data is not None:
            if np.nanmax(data) > vmax and np.nanmin(data) < vmin:
                extend_arg = 'both'
            elif np.nanmax(data) > vmax:
                extend_arg = 'max'
            elif np.nanmin(data) < vmin:
                extend_arg = 'min'
            
        norm = mpl.colors.BoundaryNorm(bounds, colorbar_cmap.N, extend=extend_arg)
    if matplot_axes is not None:
        fig = matplot_axes.figure
        ax = matplot_axes
    else:
        fig, ax = plt.subplots(figsize=figsize, layout='constrained')
    cbar_obj=fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=colorbar_cmap),cax=ax, orientation=orientation,label=cbar_label,**kwargs)
    if cbar_title:
        ax.set_title(cbar_title, fontsize='small')
    if tick_label:
        cbar_obj.set_ticklabels(tick_label, fontsize='small')
    if show_plot == True:
        plt.rcParams['figure.dpi'] = image_res
        plt.show()
    if save_to_file == True:
        plt.rcParams['savefig.dpi'] = image_res
        if output_filepath is None:
            output_filepath = f'{os.path.dirname(os.path.abspath(__file__))}/colorbar.png'
        plt.savefig(output_filepath)  
#%%
def plot_annotation(annot: list|np.ndarray|pd.Series|pd.DataFrame,
                    cmap: str|ListedColormap|LinearSegmentedColormap|None=None,
                    figsize:list=[2,2],
                    show_plot:bool = False, 
                    save_to_file: bool = False,
                    output_filepath: str | None = None, 
                    image_res: int = 600, 
                    matplot_axes = None
                    ):
    if isinstance(annot, list|np.ndarray):
        annot_plot = np.array(annot)
        annot_plot.sort()
        annot_uniq=np.unique(annot_plot)
    elif isinstance(annot, pd.Series):
        annot_plot = annot.values
        annot_uniq=annot.sort_values().unique()
    if isinstance(cmap, ListedColormap|LinearSegmentedColormap):
        print('check1')
        annot_cmap = cmap
    else:
        cmap = cmap or 'Blues'
        ncolors = 256
        annot_bound=np.concatenate(([-0.5], annot_uniq))
        annot_cmap= plt.get_cmap(cmap, len(annot_bound))
        norm = BoundaryNorm(annot_bound, len(annot_uniq))

    ax = matplot_axes if matplot_axes is not None else plt.subplots(figsize=figsize)[1]
    age_annot = ax.imshow(annot_plot.reshape(-1, 1), aspect = 'auto',cmap= annot_cmap,)
    ax.margins(x=0, y=0)
    if show_plot == True:
        plt.rcParams['figure.dpi'] = image_res
        plt.show()
    if save_to_file == True:
        plt.rcParams['savefig.dpi'] = image_res
        if output_filepath is None:
            output_filepath = f'{os.path.dirname(os.path.abspath(__file__))}/annot_bar.png'
        plt.savefig(output_filepath)  
#%%
def format_func(value, tick_number):
    if value < 0.1:
        return f"{value:.0e}"
    else:
        return f"{value:.2f}"
formatter = ticker.FuncFormatter(format_func)
#%%
def plot_basecount(alignment:np.ndarray|None,
                   figsize:list=[2,2],
                   matplot_axes=None,
                   bar_title: str|None = None,
                   xlim: list|None = None,
                   ylim: list|None = None,
                   show_plot:bool = False,
                   save_to_file: bool = False,
                   image_res: int = 600,
                   opacity:int|None=None,
                   output_filepath: str | None = None,
                   **kwargs):
    col_dict_dna={'-': 'white','A':'green','C': 'blue','T': 'red','G': 'yellow',}
    base_count=mapper.base_count(alignment=alignment)
    ax = matplot_axes if matplot_axes is not None else plt.subplots(figsize=figsize)[1]
    bottom = np.zeros(base_count.shape[1])
    for values, label in zip(base_count, ['-','A','C','T','G']):
        ax.bar(np.arange(base_count.shape[1]), values, label=label, bottom=bottom, color=col_dict_dna[label], **kwargs)
        bottom += values
        ax.margins(x=0, y=0)
    if bar_title is not None:
        ax.set_title(bar_title)
    if xlim is not None:
        ax.set_xlim(xlim[0],xlim[1])
    if ylim is not None:
        ax.set_ylim(ylim[0],ylim[1])
    if show_plot == True:
        plt.rcParams['figure.dpi'] = image_res
        plt.show()
    if save_to_file == True:
        plt.rcParams['savefig.dpi'] = image_res
        if output_filepath is None:
            output_filepath = f'{os.path.dirname(os.path.abspath(__file__))}/barplot.png'
        plt.savefig(output_filepath)  
#%%
def plot_logos(alignment:np.ndarray|None,
                figsize:list=[2,2],
                matplot_axes=None,
                bar_title: str|None = None,
                xlim: list|None = None,
                ylim: list|None = None,
                show_plot:bool = False,
                save_to_file: bool = False,
                image_res: int = 600,
                output_filepath: str | None = None,
                opacity:int|None=None,
                highlights: list|None = None,
                highlight_col: list|None = None,
                **kwargs):
    import logomaker
    base_count=mapper.base_count(alignment=alignment)
    bese_count_nogap=base_count[1:]
    counts_mat=bese_count_nogap.transpose()
    column_names = ['A', 'C', 'T', 'G']
    counts_mat_df = pd.DataFrame(counts_mat, columns=column_names)
    info_mat = logomaker.transform_matrix(counts_mat_df, from_type='counts', to_type='information')
    ax = matplot_axes if matplot_axes is not None else plt.subplots(figsize=figsize)[1]
    logos=logomaker.Logo(info_mat, ax=ax, **kwargs)
    if highlights is not None:
        for idx,highlight in enumerate(highlights):
            logos.highlight_position_range(pmin=highlight[0], pmax=highlight[1], color = highlight_col[idx])
    ax.margins(x=0, y=0)
    if bar_title is not None:
        ax.set_title(bar_title)
    if xlim is not None:
        ax.set_xlim(xlim[0],xlim[1])
    if ylim is not None:
        ax.set_ylim(ylim[0],ylim[1])
    if show_plot == True:
        plt.rcParams['figure.dpi'] = image_res
        plt.show()
    if save_to_file == True:
        plt.rcParams['savefig.dpi'] = image_res
        if output_filepath is None:
            output_filepath = f'{os.path.dirname(os.path.abspath(__file__))}/barplot.png'
        plt.savefig(output_filepath)  
#%%
_PLOT_LAYOUT = Literal['overlay','spread']
def all_overlay_plot(data: list|np.ndarray|None = None,
                     alignment: np.ndarray|None = None,
                    save_to_file: bool = False,
                    heatmap:bool = True,
                    aggregated:bool = False,
                    heatmap_annot: pd.DataFrame|None= None,
                    colorbar:bool = False,
                    logos:bool = False,
                    base_count:bool = False,
                    output_filepath: str | None = None, 
                    image_res: int = 600, 
                    alignment_col: str = 'nulc_white',
                    h_cmap: ListedColormap|str|None = None,
                    a_colset: str|None = None,
                    cbar_steps: list|None = None,
                    vlim: list|None = None,
                    figsize:list=[10,10],
                    show_plot:bool = True,
                    **kwargs):
    plt.rcParams['savefig.dpi'] = image_res
    plt.rcParams['figure.dpi'] = image_res
    fig = plt.figure(figsize=figsize)
    grid = fig.add_gridspec(nrows = 100, ncols = 100, hspace=0)
    global_kwargs = {key: value for key, value in kwargs.items() if not key.startswith(('ag_', 'hm_', 'cb_','bc_','lg_'))}
    ag_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('ag_')}
    hm_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('hm_')}
    cb_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('cb_')}
    bc_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('bc_')}
    lg_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('lg_')}
    if heatmap:
        heatmap_grid = fig.add_subplot(grid[0:50,20:70])
        if alignment is not None:
            plot_heatmap(alignment, cmap = alignment_col,matplot_axes = heatmap_grid,vmin=0, vmax=5,**hm_kwargs, **global_kwargs)
        if data is not None:
            for idx, overlay in enumerate(data):
                plot_heatmap(overlay, cmap = h_cmap[idx],matplot_axes = heatmap_grid,vmin=vlim[idx][0], vmax=vlim[idx][1],**hm_kwargs,**global_kwargs)
        heatmap_grid.tick_params(axis='both', which='major', labelsize=8)
        heatmap_grid.xaxis.set_major_locator(MaxNLocator(integer=True))
        heatmap_grid.set_yticks([])
    if aggregated:
        aggregated_grid = fig.add_subplot(grid[55:65,20:70])
        for idx, overlay in enumerate(data):
            plot_bar(overlay,color = a_colset[idx] ,alignment=alignment,mode='average', matplot_axes=aggregated_grid,**ag_kwargs,**global_kwargs)
        aggregated_grid.tick_params(axis='both', which='major', labelsize=8)
    if base_count:
        base_count_grid = fig.add_subplot(grid[75:85,20:70])
        plot_basecount(alignment=alignment, matplot_axes=base_count_grid, **bc_kwargs,**global_kwargs)
        base_count_grid.tick_params(axis='both', which='major', labelsize=8)
    if logos:
        logos_grid = fig.add_subplot(grid[65:75,20:70])
        plot_logos(alignment=alignment, matplot_axes=logos_grid, **lg_kwargs,**global_kwargs)
        logos_grid.tick_params(axis='both', which='major', labelsize=8)
    if colorbar:
        cbar_grids = []
        for idx in range(len(data)):
            cbar_grids.append(fig.add_subplot(grid[0:50, 70+7*idx:72+7*idx]))
        for idx, overlay in enumerate(data):
            plot_colorbar(cmap=h_cmap[idx], data=overlay, vmin=vlim[idx][0], vmax=vlim[idx][1], matplot_axes=cbar_grids[idx], step = cbar_steps[idx],**cb_kwargs)
            cbar_grids[idx].tick_params(axis='both', which='major', labelsize=8,)
    if heatmap_annot is not None:         
        anno_grid = fig.add_subplot(grid[0:50,18:20])
        anno_grid.tick_params(axis='both', which='major', labelsize=8)
        anno_grid.set_xticks([])
        plot_annotation(heatmap_annot, show_plot=True,matplot_axes=anno_grid)
    if show_plot == True:
        plt.rcParams['figure.dpi'] = image_res
        plt.show()
    if save_to_file == True:
        plt.rcParams['savefig.dpi'] = image_res
        if output_filepath is None:
            output_filepath = f'{os.path.dirname(os.path.abspath(__file__))}/annot_bar.png'
        plt.savefig(output_filepath)  
#%%
def heatmap_overlay_agg_spread_plot(data: list|np.ndarray|None = None,
                     alignment: np.ndarray|None = None,
                    save_to_file: bool = False,
                    heatmap:bool = True,
                    aggregated:bool = False,
                    heatmap_annot: pd.DataFrame|None= None,
                    colorbar:bool = False,
                    logos:bool = False,
                    base_count:bool = False,
                    output_filepath: str | None = None, 
                    image_res: int = 600, 
                    alignment_col: str = 'nulc_white',
                    h_cmap: ListedColormap|str|None = None,
                    a_colset: str|None = None,
                    cbar_steps: list|None = None,
                    agg_titles: list|None = None,
                    vlim: list|None = None,
                    figsize:list=[10,10],
                    show_plot:bool = True,
                    **kwargs):
    plt.rcParams['savefig.dpi'] = image_res
    plt.rcParams['figure.dpi'] = image_res
    fig = plt.figure(figsize=figsize)
    grid = fig.add_gridspec(nrows = 100, ncols = 100, hspace=0)
    global_kwargs = {key: value for key, value in kwargs.items() if not key.startswith(('ag_', 'hm_', 'cb_','bc_','lg_'))}
    ag_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('ag_')}
    hm_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('hm_')}
    cb_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('cb_')}
    bc_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('bc_')}
    lg_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('lg_')}
    if heatmap:
        heatmap_grid = fig.add_subplot(grid[0:50,20:70])
        if alignment is not None:
            plot_heatmap(alignment, cmap = alignment_col,matplot_axes = heatmap_grid,vmin=0, vmax=5,**hm_kwargs, **global_kwargs)
        if data is not None:
            for idx, overlay in enumerate(data):
                plot_heatmap(overlay, cmap = h_cmap[idx],matplot_axes = heatmap_grid,vmin=vlim[idx][0], vmax=vlim[idx][1],**hm_kwargs,**global_kwargs)
        heatmap_grid.tick_params(axis='both', which='major', labelsize=8)
        heatmap_grid.xaxis.set_major_locator(MaxNLocator(integer=True))
        heatmap_grid.set_yticks([])
    if aggregated:
        aggregated_grids = []
        for idx in range(len(data)):
            aggregated_grids.append(fig.add_subplot(grid[50+10*idx:60+10*idx, 20:70]))
        for idx, overlay in enumerate(data):
            plot_bar(overlay,color = a_colset[idx] ,alignment=alignment,mode='average', matplot_axes=aggregated_grids[idx],**ag_kwargs,**global_kwargs)
            aggregated_grids[idx].tick_params(axis='both', which='major', labelsize=8)
        if agg_titles is not None:
            for idx, agg_title in enumerate(agg_titles):
                aggregated_grids[idx].set_title(agg_title, fontsize='small', x=1.05, y=0.4)
    if base_count:
        base_count_grid = fig.add_subplot(grid[50+10*len(data):60+10*len(data),20:70])
        plot_basecount(alignment=alignment, matplot_axes=base_count_grid, **bc_kwargs,**global_kwargs)
        base_count_grid.tick_params(axis='both', which='major', labelsize=8)
    if logos:
        logos_grid = fig.add_subplot(grid[50+10*(len(data)+1):60+10*(len(data)+1),20:70])
        plot_logos(alignment=alignment, matplot_axes=logos_grid, **lg_kwargs,**global_kwargs)
        logos_grid.tick_params(axis='both', which='major', labelsize=8)
    if colorbar:
        cbar_grids = []
        for idx in range(len(data)):
            cbar_grids.append(fig.add_subplot(grid[0:50, 70+7*idx:72+7*idx]))
        for idx, overlay in enumerate(data):
            plot_colorbar(cmap=h_cmap[idx], data=overlay, vmin=vlim[idx][0], vmax=vlim[idx][1], matplot_axes=cbar_grids[idx], step = cbar_steps[idx],**cb_kwargs)
            cbar_grids[idx].tick_params(axis='both', which='major', labelsize=8,)
    if heatmap_annot is not None:         
        anno_grid = fig.add_subplot(grid[0:50,18:20])
        anno_grid.tick_params(axis='both', which='major', labelsize=8)
        anno_grid.set_xticks([])
        plot_annotation(heatmap_annot, show_plot=True,matplot_axes=anno_grid)
    if show_plot == True:
        plt.rcParams['figure.dpi'] = image_res
        plt.show()
    if save_to_file == True:
        plt.rcParams['savefig.dpi'] = image_res
        if output_filepath is None:
            output_filepath = f'{os.path.dirname(os.path.abspath(__file__))}/annot_bar.png'
        plt.savefig(output_filepath)  

# %%
_SPREAD_MODE = Literal['horizontal','vertical']
def heatmap_spread_plot(data: list|np.ndarray|None = None,
                     alignment: np.ndarray|None = None,
                    save_to_file: bool = False,
                    heatmap:bool = True,
                    aggregated:bool = False,
                    heatmap_annot: pd.DataFrame|None= None,
                    colorbar:bool = False,
                    logos:bool = False,
                    base_count:bool = False,
                    output_filepath: str | None = None, 
                    image_res: int = 600, 
                    alignment_col: str = 'nulc_white',
                    h_cmap: ListedColormap|str|None = None,
                    a_colset: str|None = None,
                    cbar_steps: list|None = None,
                    vlim: list|None = None,
                    figsize:list=[50,50],
                    show_plot:bool = True,
                    spread_mode:_SPREAD_MODE = 'horizontal',
                    **kwargs):
    plt.rcParams['savefig.dpi'] = image_res
    plt.rcParams['figure.dpi'] = image_res
    fig = plt.figure(figsize=figsize)
    grid = fig.add_gridspec(nrows = 500, ncols = 500, hspace=0)
    global_kwargs = {key: value for key, value in kwargs.items() if not key.startswith(('ag_', 'hm_', 'cb_','bc_','lg_'))}
    ag_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('ag_')}
    hm_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('hm_')}
    cb_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('cb_')}
    bc_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('bc_')}
    lg_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('lg_')}
    total_data_count = len(data)
    if alignment is not None:
        total_data_count += 1
    heatmap_grids = []
    for idx in range(total_data_count):
        if spread_mode == 'horizontal':
            heatmap_grids.append(fig.add_subplot(grid[0:50,20+60*idx:70+60*idx]))
        elif spread_mode == 'vertical':
            heatmap_grids.append(fig.add_subplot(grid[0+65*idx:50+65*idx,20:70]))
    if alignment is not None:
        plot_heatmap(alignment, cmap = alignment_col,matplot_axes = heatmap_grids[0],vmin=0, vmax=5,**hm_kwargs, **global_kwargs)
        heatmap_grids[0].tick_params(axis='both', which='major', labelsize=8)
        heatmap_grids[0].xaxis.set_major_locator(MaxNLocator(integer=True))
        heatmap_grids[0].set_yticks([])
    if data is not None:
        for idx, overlay in enumerate(data):
            if alignment is not None:
                grid_idx = idx+1
            else:
                grid_idx = idx
            plot_heatmap(overlay, cmap = h_cmap[idx],matplot_axes = heatmap_grids[grid_idx],vmin=vlim[idx][0], vmax=vlim[idx][1],**hm_kwargs,**global_kwargs)
            heatmap_grids[grid_idx].tick_params(axis='both', which='major', labelsize=8)
            heatmap_grids[grid_idx].xaxis.set_major_locator(MaxNLocator(integer=True))
            heatmap_grids[grid_idx].set_yticks([])
    if aggregated:
        aggregated_grids = []
        for idx, overlay in enumerate(data):
            if alignment is not None:
                grid_idx = idx+1
            else:
                grid_idx = idx
            if spread_mode == 'horizontal':
                aggregated_grid = fig.add_subplot(grid[50:60,20+60*grid_idx:70+60*grid_idx])
            elif spread_mode == 'vertical':
                aggregated_grid = fig.add_subplot(grid[50+65*grid_idx:60+65*grid_idx,20:70])
            plot_bar(overlay,color = a_colset[idx] ,alignment=alignment,mode='average', matplot_axes=aggregated_grid,**ag_kwargs,**global_kwargs)
            aggregated_grid.tick_params(axis='both', which='major', labelsize=8)
    if base_count:
        base_count_grid = fig.add_subplot(grid[50:60,20:70])
        plot_basecount(alignment=alignment, matplot_axes=base_count_grid, **bc_kwargs,**global_kwargs)
        base_count_grid.tick_params(axis='both', which='major', labelsize=8)
    if logos:
        logos_grid = fig.add_subplot(grid[50:60,20:70])
        plot_logos(alignment=alignment, matplot_axes=logos_grid, **lg_kwargs,**global_kwargs)
        logos_grid.tick_params(axis='both', which='major', labelsize=8)
    if colorbar:
        cbar_grids = []
        for idx in range(len(data)):
            if alignment is not None:
                grid_idx = idx+1
            else:
                grid_idx = idx
            if spread_mode == 'horizontal':
                cbar_grids.append(fig.add_subplot(grid[0:50,70+60*grid_idx:72+60*grid_idx]))
            elif spread_mode == 'vertical':
                cbar_grids.append(fig.add_subplot(grid[0+65*grid_idx:50+65*grid_idx,70:72]))
        for idx, overlay in enumerate(data):
            plot_colorbar(cmap=h_cmap[idx], data=overlay, vmin=vlim[idx][0], vmax=vlim[idx][1], matplot_axes=cbar_grids[idx], step = cbar_steps[idx],**cb_kwargs)
            cbar_grids[idx].tick_params(axis='both', which='major', labelsize=8,)
    if heatmap_annot is not None:         
        for idx in range(total_data_count):
            if spread_mode == 'horizontal':
                anno_grid = fig.add_subplot(grid[0:50,18+60*idx:20+60*idx])
                anno_grid.set_xticks([])
                anno_grid.tick_params(axis='both', which='major', labelsize=8)
            elif spread_mode == 'vertical':
                anno_grid = fig.add_subplot(grid[0+65*idx:50+65*idx,18:20])
                anno_grid.set_xticks([])
                anno_grid.tick_params(axis='both', which='major', labelsize=8)
            plot_annotation(heatmap_annot, show_plot=True,matplot_axes=anno_grid)
            
            
    if show_plot == True:
        plt.rcParams['figure.dpi'] = image_res
        plt.show()
    if save_to_file == True:
        plt.rcParams['savefig.dpi'] = image_res
        if output_filepath is None:
            output_filepath = f'{os.path.dirname(os.path.abspath(__file__))}/annot_bar.png'
        plt.savefig(output_filepath)