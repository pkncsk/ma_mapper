#%%
from configparser import Interpolation
import sys
from tabnanny import verbose
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import mapper
#%%
subfamily = ['THE1C']
alignment_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/alignment/'+subfamily[0]+'.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered, reference_table='/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/combined_age_div/combined_te_age_div.txt')
# %%
coord_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/coord_internal_id/'+subfamily[0]+'.txt'
#%%
bam_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/KZFP-bam_hg38/znf267.sorted.bam'
bam_forward=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_forward')
bam_reverse=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_reverse')
bam_min=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_min')
#%%
bigwig_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_bigwig/241-mammalian-2020v2.bigWig'
phylop=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig')
#%%
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/zoonomia447/hg38.phyloP447way.bw'
phylop_447=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig')
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/zoonomia447/hg38.phyloP447wayLRT.bw'
phylop_LRT=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig')
#%%
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/AP-1(bZIP).bed'
ap1=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed')
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/NFkB-p65.bed'
nfkb=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed')
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/BMAL1.bed'
bmal=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed')
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/CLOCK.bed'
clock=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed')
#%%

# %%
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
if sys.version_info >= (3, 8, 0):
    from typing import Literal, Tuple, List
else:
    from typing_extensions import Literal, Tuple, List
#
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
             plot_title: str|None = None,
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
import pandas as pd
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
        annot_cmap = cmap
    else:
        cmap = cmap or 'Blues'
        ncolors = 256
        annot_bound=np.concatenate(([-0.5], annot_uniq))
        annot_cmap= plt.get_cmap(cmap, len(annot_bound))
        norm = BoundaryNorm(annot_bound, len(annot_uniq))

    ax = matplot_axes if matplot_axes is not None else plt.subplots(figsize=figsize)[1]
    age_annot = ax.imshow(annot_plot.reshape(-1, 1), aspect = 'auto',cmap= annot_cmap, interpolation='nearest', norm=norm)
    if show_plot == True:
        plt.rcParams['figure.dpi'] = image_res
        plt.show()
    if save_to_file == True:
        plt.rcParams['savefig.dpi'] = image_res
        if output_filepath is None:
            output_filepath = f'{os.path.dirname(os.path.abspath(__file__))}/annot_bar.png'
        plt.savefig(output_filepath)  
#%%
#logos imitation 
alignment_transposed = alignment_filtered.transpose()
output_array = []
for idx, ref_count in enumerate(alignment_transposed):
    (unique, counts) = np.unique(ref_count, return_counts=True)
    frequencies = dict(zip(unique, counts))
    base_count = []
    for base in [0,1,2,3,4]:
        nucl_count = frequencies.get(base,0)
        base_count.append(nucl_count)
    output_array.append(base_count)
base_count=np.array(output_array).transpose()
#%%
col_dict_dna={
    '-': 'white',
    'A':'green',
    'C': 'blue',
    'T': 'red',
    'G': 'yellow',
}
#%%
import logomaker
bese_count_nogap=base_count[1:]
#logomaker need A  C  G  T
#our format is A C T G, need to switch data entries
#bese_count_nogap[[2,3]] = bese_count_nogap[[3,2]]
#transpose
counts_mat=bese_count_nogap.transpose()
column_names = ['A', 'C', 'T', 'G']
counts_mat_df = pd.DataFrame(counts_mat, columns=column_names)
info_mat = logomaker.transform_matrix(counts_mat_df, from_type='counts', to_type='information')
#logomaker.Logo(info_mat)
#plot_annotation(metadata_age.te_age, show_plot=True,)
# %%
from matplotlib.ticker import MaxNLocator
xlim = None
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['figure.dpi'] = 600
fig = plt.figure(figsize=(10,10))
grid = fig.add_gridspec(nrows = 100, ncols = 100, hspace=0)
upp = fig.add_subplot(grid[0:50,20:70])
dow = fig.add_subplot(grid[55:65,20:70])
dow2 = fig.add_subplot(grid[65:75,20:70])
dow3 = fig.add_subplot(grid[75:85,20:70])
cbar = fig.add_subplot(grid[0:50, 70:72])
cbar2 = fig.add_subplot(grid[0:50, 77:79])
cbar3 = fig.add_subplot(grid[0:50, 84:86])
cbar4 = fig.add_subplot(grid[0:50, 91:93])
cbar5 = fig.add_subplot(grid[50:70, 71:73])
anno = fig.add_subplot(grid[0:50,18:20])
upp.tick_params(axis='both', which='major', labelsize=8)
dow.tick_params(axis='both', which='major', labelsize=8)
cbar.tick_params(axis='both', which='major', labelsize=8)
cbar2.tick_params(axis='both', which='major', labelsize=8)
cbar3.tick_params(axis='both', which='major', labelsize=8)
cbar4.tick_params(axis='both', which='major', labelsize=8)
cbar5.tick_params(axis='both', which='major', labelsize=8)
anno.tick_params(axis='both', which='major', labelsize=8)
from ma_mapper import custom_cmap
plot_heatmap(phylop_447, cmap=custom_cmap.vlag_r_mpl,vmax=0.5,vmin=-0.5, matplot_axes = upp, zero_col_centering=False, xlim=xlim)
upp.xaxis.set_major_locator(MaxNLocator(integer=True))
upp.set_yticks([])
plot_heatmap(ap1, cmap = 'Greens',matplot_axes = upp,vmin=0, vmax=1, opacity=0.5, xlim=xlim)
plot_heatmap(bam_forward, cmap = 'Reds',matplot_axes = upp,vmin=0, vmax=0.001, xlim=xlim)
plot_heatmap(bam_reverse, cmap = 'Blues',matplot_axes = upp,vmin=0, vmax=0.001, xlim=xlim)
#upp.set_xticks([])
anno.set_xticks([])
bottom = np.zeros(base_count.shape[1])
for values, label in zip(base_count, ['-','A','C','T','G']):
    dow3.bar(np.arange(base_count.shape[1]), values, label=label, bottom=bottom, color=col_dict_dna[label])
    bottom += values
#dow3.set_xlim(xlim[0],xlim[1])
logos=logomaker.Logo(info_mat, ax=dow2)
logos.highlight_position_range(pmin=293, pmax=299, color='grey')
#dow2.set_xlim(xlim[0],xlim[1])
plot_bar(phylop_447, alignment=alignment_filtered,mode='average', matplot_axes=dow, xlim=xlim)
plot_bar(ap1, alignment=alignment_filtered,mode='average', matplot_axes=dow, color='green', ylim=[-0.5,0.5], opacity=0.5, xlim=xlim)
plot_bar(bam_forward, alignment=alignment_filtered,mode='average', matplot_axes=dow, color='red', ylim=[-0.5,0.5], opacity=0.5, xlim=xlim)
plot_bar(bam_reverse, alignment=alignment_filtered,mode='average', matplot_axes=dow, color='blue', ylim=[-0.5,0.5], opacity=0.5, xlim=xlim)
plot_colorbar(cmap=custom_cmap.vlag_r_mpl,data=phylop_447, vmax=0.5,vmin=-0.5,matplot_axes=cbar, step = 0.1)
plot_colorbar(cmap='Greens',data=ap1,matplot_axes=cbar2, step = 0.1,vmin=0,vmax=1)
import matplotlib.ticker as ticker
def format_func(value, tick_number):
    if value < 0.1:
        return f"{value:.0e}"
    else:
        return f"{value:.2f}"
formatter = ticker.FuncFormatter(format_func)

plot_colorbar(cmap='Reds',data=bam_forward,matplot_axes=cbar3, step = 0.0001,vmin=0,vmax=0.001)
plot_colorbar(cmap='Blues',data=bam_forward,matplot_axes=cbar4, step = 0.0001,vmin=0,vmax=0.001)
plot_colorbar(cmap='Blues',data=bam_forward,matplot_axes=cbar4, step = 0.0001,vmin=0,vmax=0.001)
plot_colorbar(cmap='dna_jaspar',data=alignment_filtered,matplot_axes=cbar5, tick_label= ['gap','A','C','T','G'])
cbar3.yaxis.set_major_formatter(formatter)
cbar4.yaxis.set_major_formatter(formatter)
upp.margins(x=0, y=0)
dow3.margins(x=0, y=0)
dow3.xaxis.set_major_locator(MaxNLocator(integer=True))
plot_annotation(metadata_age.te_age, show_plot=True,matplot_axes=anno)
plt.show()
#%%
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import plots
#%%
import importlib
importlib.reload(plots)
plots.all_overlay_plot(alignment=alignment_filtered, alignment_col='nulc_white',data=[phylop_447, ap1],h_cmap=[custom_cmap.vlag_r_mpl, 'Greens'],vlim=[[-0.5,0.5],[0,1]], aggregated=True, a_colset = ['grey','green'], colorbar=True, cbar_steps=[0.1,0.1], heatmap_annot=metadata_age.te_age)
#%%
plot_heatmap(alignment_filtered, cmap = 'dna',vmin=0, vmax=6, opacity=0.5, Interpolation='nearest',)
#%%
from matplotlib.ticker import MaxNLocator
xlim = None
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['figure.dpi'] = 600
fig = plt.figure(figsize=(10,10))
grid = fig.add_gridspec(nrows = 100, ncols = 100, hspace=0)
upp = fig.add_subplot(grid[0:50,20:70])
dow = fig.add_subplot(grid[55:65,20:70])

cbar = fig.add_subplot(grid[0:50, 70:72])

anno = fig.add_subplot(grid[0:50,18:20])
upp.tick_params(axis='both', which='major', labelsize=8)
dow.tick_params(axis='both', which='major', labelsize=8)
cbar.tick_params(axis='both', which='major', labelsize=8)
anno.tick_params(axis='both', which='major', labelsize=8)
from ma_mapper import custom_cmap
plot_heatmap(phylop_LRT, cmap=custom_cmap.vlag_r_mpl,vmax=0.5,vmin=-0.5, matplot_axes = upp, zero_col_centering=False, xlim=xlim)
upp.xaxis.set_major_locator(MaxNLocator(integer=True))
upp.set_yticks([])
#upp.set_xticks([])
anno.set_xticks([])
plot_bar(phylop_LRT, alignment=alignment_filtered,mode='average', matplot_axes=dow, xlim=xlim)
plot_colorbar(cmap=custom_cmap.vlag_r_mpl,data=phylop_LRT, vmax=0.5,vmin=-0.5,matplot_axes=cbar, step = 0.1)

upp.margins(x=0, y=0)

plot_annotation(metadata_age.te_age, show_plot=True,matplot_axes=anno)
plt.show()

#%%
import matplotlib as mpl
fig, ax = plt.subplots(figsize=(10,1), layout='constrained')
cmap = custom_cmap.vlag_r_mpl
bounds =np.arange(-0.5,0.6,0.1)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             cax=ax, orientation='horizontal',
             label="Discrete intervals with extend='both' keyword")
#%%
fig, ax = plt.subplots(figsize=(10,1), layout='constrained')

#%%
from ma_mapper._unused import plot
plot.overlay_plot([alignment_filtered,phylop,bam_forward, bam_reverse,ap1,nfkb,bmal,clock], metadata= metadata_age, nucleotide_color='white',data_vmin=[-0.5,0,0,0,0,0,0],data_vmax = [0.5,0.1,0.1,1,1,1,1],data_cmap=['RdBu','Purples','Oranges','Blues','Reds','Greens','copper'], plot_title='ChIPexo_znf276_tfbs', show_data_legend=True,data_label=['phyloP','ChIP_forward','ChIP_reverse','AP1','nfkB-p65','BMAL','CLOCK'])
#%%
#%%
from ma_mapper import mapper

phylop_mean=mapper.normalise(alignment=alignment_filtered,mapped_data=phylop)
fig, ax = plt.subplots(figsize=(10,2))
ax.bar(range(len(phylop_mean)), phylop_mean, color = 'grey')
ax.margins(x=0, y=0)

#%%
fig, ax = plt.subplots(figsize=(2,2))
nucleotide_labels = ['gap', 'A', 'C', 'T', 'G']
nucleotide_color_list = ['grey','green','yellow','red','blue']
align_cmap = ListedColormap(nucleotide_color_list)

#align_cmap =ListedColormap(['grey','white','white','white','white'])
ax.set_alpha(0.5)
heatmap = ax.imshow(alignment_filtered, aspect = 'auto',cmap= align_cmap, interpolation='nearest', vmin=0,vmax=5)

#%%
fig, ax = plt.subplots(figsize=(5,5))
heatmap = ax.imshow(phylop,aspect='auto',cmap= custom_cmap.vlag_r_mpl, vmax=0.5, vmin=-0.5)

heatmap = ax.imshow(bam_min,aspect='auto',cmap= 'Greens', vmax=0.5, vmin=-0.5)
# %%
fig, ax = plt.subplots(figsize=(2,2))
ax.set_alpha(0.1)
heatmap = ax.imshow(alignment_filtered, aspect = 'auto',cmap= align_cmap, interpolation='nearest', vmin=0,vmax=5)
ax.set_alpha(0.9)
heatmap = ax.imshow(phylop,aspect='auto',cmap= overlay_cmap_alpha, interpolation='nearest', vmax=0.5, vmin=-0.5)
# %%
from ma_mapper import custom_cmap
plot_heatmap(phylop, cmap=custom_cmap.vlag_r_mpl, image_res=600, vmax=0.5,vmin=-0.5,figsize=[4,4], zero_col_centering =True)
#%%
plot_bar(phylop, alignment=alignment_filtered,mode='average', image_res=600,figsize=[10,2])
#