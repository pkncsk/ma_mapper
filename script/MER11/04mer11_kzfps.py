#%%
import sys
from turtle import width
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
from ma_mapper import mapper
#%%
subfamily = 'MER11'
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file, col_threshold = 0.5, row_threshold = 0.5)
age_table_list = [f'{config.te_age_folder}/MER11A.txt',
             f'{config.te_age_folder}/MER11B.txt',
             f'{config.te_age_folder}/MER11C.txt']
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered, age_table=age_table_list)
metadata_age['len'] = metadata_age.end.astype(int) - metadata_age.start.astype(int)
# %%
coord_file = f'{config.coord_internal_id_folder}/{subfamily}.txt'
#%%
bam_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/KZFP-bam_hg38/znf808.sorted.bam'

znf808_forward=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_forward', pf_col_threshold = 0.5, pf_row_threshold = 0.5)
znf808_reverse=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_reverse', pf_col_threshold = 0.5, pf_row_threshold = 0.5)
#%%
bam_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/KZFP-bam_hg38/znf525.sorted.bam'
znf525_forward=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_forward', pf_col_threshold = 0.5, pf_row_threshold = 0.5)
znf525_reverse=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_reverse', pf_col_threshold = 0.5, pf_row_threshold = 0.5)
#%%
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/zoonomia447/hg38.phyloP447way.bw'
phylop_447=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', pf_col_threshold = 0.5, pf_row_threshold = 0.5)
#%%
from ma_mapper import plots
plots.all_overlay_plot(data = [znf808_forward,znf808_reverse,znf525_forward,znf525_reverse], alignment= alignment_filtered,heatmap_annot=metadata_age.te_age, h_cmap=['Blues','Reds','Greens','Oranges'], vlim = [[0,0.1],[0,0.1],[0,0.01],[0,0.01]], opacity = 0.5, aggregated=True,a_colset = ['blue','red','green','orange'])
#%%
import numpy as np
znf808_minimum = np.minimum(np.nanmean(znf808_forward, axis = 0), np.nanmean(znf808_reverse, axis = 0))
znf525_minimum = np.minimum(np.nanmean(znf525_forward, axis = 0), np.nanmean(znf525_reverse, axis = 0))
# %%
import importlib
importlib.reload(plots)
plots.heatmap_overlay_agg_spread_plot(data=[znf808_minimum,znf525_minimum], heatmap=False, opacity=0.5, a_colset = ['purple','brown'], aggregated=True, mode='1D')
# %%
from scipy.signal import find_peaks
peaks, _ = find_peaks(znf808_minimum, width = 8)
highest_peak_index = peaks[np.argmax(znf808_minimum[peaks])]
forward_nonzero_indices = np.where(znf808_forward[:, highest_peak_index] != 0)[0]
reverse_nonzero_indices = np.where(znf808_reverse[:, highest_peak_index] != 0)[0]
znf808_binding_indices=list(set(forward_nonzero_indices) | set(reverse_nonzero_indices))
znf808_nonbinding_indices=list(set(np.arange(len(znf808_forward))) - set(znf808_binding_indices))
#%% big figure ZNF808
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from ma_mapper import plots
import importlib
import matplotlib as mpl
import matplotlib.colors as mcolors
importlib.reload(plots)
xlim = None
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['figure.dpi'] = 600
fig = plt.figure(figsize=(10,10))
#cbar2 = fig.add_subplot(grid[0:50, 77:79])
#cbar2.tick_params(axis='both', which='major', labelsize=8)
from ma_mapper import custom_cmap

grid = fig.add_gridspec(nrows = 100, ncols = 100, hspace=0)
total=len(phylop_447)
nonzero_grid=round(len(znf808_binding_indices)/total * 50)
heatmap_g1 = fig.add_subplot(grid[0:nonzero_grid,20:70])
###########################################################
phylop_447_bind=phylop_447[znf808_binding_indices]
alignemnt_bind=alignment_filtered[znf808_binding_indices]
phylop_447_bind[alignemnt_bind == 0] = np.nan
phylop_447_nonbind=phylop_447[znf808_nonbinding_indices]
alignemnt_nonbind=alignment_filtered[znf808_nonbinding_indices]
phylop_447_nonbind[alignemnt_nonbind == 0] = np.nan
##########################################################
plots.plot_heatmap(phylop_447_bind, cmap=custom_cmap.vlag_r_mpl,vmax=0.5,vmin=-0.5, matplot_axes = heatmap_g1)
heatmap_g1.xaxis.set_major_locator(MaxNLocator(integer=True))
heatmap_g1.set_yticks([])
heatmap_g1.set_title('phyloP_447 MER11 grouped by ZNF808 binding')
heatmap_g1.tick_params(axis='both', which='major', labelsize=8)
plots.plot_heatmap(znf808_forward[znf808_binding_indices], cmap = 'Blues',matplot_axes = heatmap_g1,vmin=0, vmax=1, opacity=0.5,)
plots.plot_heatmap(znf808_reverse[znf808_binding_indices], cmap = 'Reds',matplot_axes = heatmap_g1,vmin=0, vmax=1, opacity=0.5,)
plots.plot_heatmap(znf525_forward[znf808_binding_indices], cmap = 'Greens',matplot_axes = heatmap_g1,vmin=0, vmax=0.1, opacity=0.5,)
plots.plot_heatmap(znf525_reverse[znf808_binding_indices], cmap = 'Oranges',matplot_axes = heatmap_g1,vmin=0, vmax=0.1, opacity=0.5,)

heatmap_g2 = fig.add_subplot(grid[nonzero_grid:50,20:70])
heatmap_g2.tick_params(axis='both', which='major', labelsize=8)
plots.plot_heatmap(phylop_447_nonbind, cmap=custom_cmap.vlag_r_mpl,vmax=0.5,vmin=-0.5, matplot_axes = heatmap_g2)
heatmap_g2.xaxis.set_major_locator(MaxNLocator(integer=True))
heatmap_g2.set_yticks([])
plots.plot_heatmap(znf808_forward[znf808_nonbinding_indices], cmap = 'Blues',matplot_axes = heatmap_g2,vmin=0, vmax=1, opacity=0.5,)
plots.plot_heatmap(znf808_reverse[znf808_nonbinding_indices], cmap = 'Reds',matplot_axes = heatmap_g2,vmin=0, vmax=1, opacity=0.5,)
plots.plot_heatmap(znf525_forward[znf808_nonbinding_indices], cmap = 'Greens',matplot_axes = heatmap_g2,vmin=0, vmax=0.1, opacity=0.5,)
plots.plot_heatmap(znf525_reverse[znf808_nonbinding_indices], cmap = 'Oranges',matplot_axes = heatmap_g2,vmin=0, vmax=0.1, opacity=0.5,)

bar = fig.add_subplot(grid[50:60,20:70])
bar.tick_params(axis='both', which='major', labelsize=8)
plots.plot_bar(phylop_447_bind, alignment=alignemnt_bind,mode='average', matplot_axes=bar,color='grey',opacity=0.5)
plots.plot_bar(phylop_447_nonbind, alignment=alignemnt_nonbind,mode='average', matplot_axes=bar,color ='grey',opacity=0.5)

bar2 = fig.add_subplot(grid[60:70,20:70])
bar2.tick_params(axis='both', which='major', labelsize=8)
plots.plot_bar(np.nanmedian(phylop_447_bind,axis=0)-np.nanmedian(phylop_447_nonbind,axis=0), alignment=alignment_filtered,mode='1D', matplot_axes=bar2)

bar3 = fig.add_subplot(grid[70:80,20:70])
bar3.tick_params(axis='both', which='major', labelsize=8)
plots.plot_bar(np.log(np.nanmedian(phylop_447_bind,axis=0)/np.nanmedian(phylop_447_nonbind,axis=0)), alignment=alignment_filtered,mode='1D', matplot_axes=bar3)

cbar = fig.add_subplot(grid[0:50, 70:72])
cbar.tick_params(axis='both', which='major', labelsize=8)
plots.plot_colorbar(cmap=custom_cmap.vlag_r_mpl,data=phylop_447, vmax=0.5,vmin=-0.5,matplot_axes=cbar, step = 0.1)
#plots.plot_colorbar(cmap='Purples',data=ap1,matplot_axes=cbar2, step = 0.1,vmin=0,vmax=1)
import matplotlib.ticker as ticker
def format_func(value, tick_number):
    if value < 0.1:
        return f"{value:.0e}"
    else:
        return f"{value:.2f}"
formatter = ticker.FuncFormatter(format_func)
heatmap_g1.margins(x=0, y=0)
heatmap_g2.margins(x=0, y=0)
######################################################age annotation

anno = fig.add_subplot(grid[0:50,18:20])
anno.tick_params(axis='both', which='major', labelsize=8)
anno.set_xticks([])
anno.set_yticks([])
cbar_anno =fig.add_subplot(grid[0:50,8:10])
cbar_anno.tick_params(axis='both', which='major', labelsize=8)
anno_sort = metadata_age.iloc[np.concatenate((znf808_binding_indices,znf808_nonbinding_indices))].te_age.fillna(0)
annot_plot = anno_sort.values
annot_uniq=anno_sort.sort_values().unique()

# Extract colors from the 'Blues' colormap
blues_cmap = plt.cm.Blues
num_colors = len(annot_uniq)
colors = [blues_cmap(i / (num_colors - 1)) for i in range(num_colors)]
# Create a colormap with a distinct color for each unique value
annot_cmap = mcolors.ListedColormap(colors)
# Add an extra boundary to center the ticks
boundaries = np.append(annot_uniq, annot_uniq[-1] + (annot_uniq[-1] - annot_uniq[-2]))
norm = mcolors.BoundaryNorm(boundaries, annot_cmap.N)

# Calculate tick positions to center them in the squares
tick_positions = (boundaries[:-1] + boundaries[1:]) / 2
# Create the plot
age_annot= anno.imshow(annot_plot.reshape(-1, 1),aspect = 'auto', cmap=annot_cmap, norm=norm)
#anno.set_title('age', fontsize='small')
# Create a colorbar with ticks at the calculated positions
anno_cbar=fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=annot_cmap),
             cax=cbar_anno, orientation='vertical', ticks= tick_positions)
anno_cbar.ax.set_yticklabels([str(val) for val in annot_uniq])
cbar_anno.set_title('(MYA)', fontsize='small')
cbar_anno.yaxis.set_ticks_position('left')
########################################################subgroup annotation
from matplotlib.colors import ListedColormap
anno2 = fig.add_subplot(grid[0:50,16:18])
anno2.tick_params(axis='both', which='major', labelsize=8)
anno2.set_xticks([])
cbar_anno2 =fig.add_subplot(grid[0:50,0:2])
cbar_anno2.tick_params(axis='both', which='major', labelsize=8)
metadata_age['subgroup'] = metadata_age['meta_id'].str.split('_').str[0]
subgroups = np.unique(metadata_age['subgroup'] )
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_annot=metadata_age['subgroup'].map(numerical_subgroup)[np.concatenate((znf808_binding_indices,znf808_nonbinding_indices))]

# Define a colormap (you can customize this)
colors = ['red', 'yellow', 'blue']  # Adjust these colors as needed
subgroup_cmap = ListedColormap(colors)
# Calculate boundaries for the categories
boundaries = np.arange(len(subgroups) + 1) - 0.5
s_norm = mcolors.BoundaryNorm(boundaries, subgroup_cmap.N)

age_annot= anno2.imshow(subgroup_annot.values.reshape(-1, 1),aspect = 'auto', cmap=subgroup_cmap)
#anno2.set_title('subgroup', fontsize='small')
anno_cbar2=fig.colorbar(mpl.cm.ScalarMappable(norm=s_norm, cmap=subgroup_cmap),
             cax=cbar_anno2, orientation='vertical', ticks= np.arange(len(subgroups)))
anno_cbar2.ax.set_yticklabels(['MER11A', 'MER11B', 'MER11C'])
cbar_anno2.set_title('subgroup', fontsize='small')
cbar_anno2.yaxis.set_ticks_position('left')
########################################################
import scipy
bar4 = fig.add_subplot(grid[80:90,20:70])
bar4.tick_params(axis='both', which='major', labelsize=8)
stat_v, p_value = scipy.stats.mannwhitneyu(phylop_447_bind,phylop_447_nonbind, axis =0,nan_policy='omit', alternative='greater')
plots.plot_bar(-np.log(p_value), alignment=alignment_filtered,mode='1D', matplot_axes=bar4)
#########################################################
bar5 = fig.add_subplot(grid[90:95,20:70])
bar5.tick_params(axis='both', which='major', labelsize=8)
def min_max_scale(lst):
    if len(lst) == 1:
        return [0.0]  # or [0.5], depending on the desired output
    min_val = min(lst)
    max_val = max(lst)
    if min_val == max_val:
        return [0.0 for _ in lst]  # or [0.5 for _ in lst]
    return [(x - min_val) / (max_val - min_val) for x in lst]
plots.plot_bar(min_max_scale(znf808_minimum), alignment=alignment_filtered,mode='1D', matplot_axes=bar5,color='purple',opacity=0.5)
plots.plot_bar(min_max_scale(znf525_minimum), alignment=alignment_filtered,mode='1D', matplot_axes=bar5,color='brown',opacity=0.5)
cbar.set_title('phyloP', fontsize='small')
#cbar2.set_title('AP1', fontsize='small')
bar.set_title('phyloP', fontsize='small', x=1.06, y=0.4)
bar2.set_title('diff', fontsize='small', x=1.03, y=0.4)
bar3.set_title('log ratio', fontsize='small', x=1.07, y=0.4)
bar4.set_title('-log p-value\nranked test\nH1=greater', fontsize='small', x=1.12, y=0.4)
bar5.set_title('KZFP binding', fontsize='small', x=1.1, y=0.4)
plt.show()
#%%
#peaks, _ = find_peaks(znf525_minimum, width=5)
peak = [732,736,747,754,779,795] #manually curated
#highest_peak_index = peaks[np.argmax(znf525_minimum[peaks])]
forward_nonzero_indices = np.where(znf525_forward[:, peaks] != 0)[0]
reverse_nonzero_indices = np.where(znf525_reverse[:, peaks] != 0)[0]
znf525_binding_indices=list(set(forward_nonzero_indices) | set(reverse_nonzero_indices))
znf525_nonbinding_indices=list(set(np.arange(len(znf525_forward))) - set(znf525_binding_indices))
#%% big figure ZNF528
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from ma_mapper import plots
import importlib
import matplotlib as mpl
import matplotlib.colors as mcolors
importlib.reload(plots)
xlim = None
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['figure.dpi'] = 600
fig = plt.figure(figsize=(10,10))
#cbar2 = fig.add_subplot(grid[0:50, 77:79])
#cbar2.tick_params(axis='both', which='major', labelsize=8)
from ma_mapper import custom_cmap
grid = fig.add_gridspec(nrows = 100, ncols = 100, hspace=0)
total=len(phylop_447)
nonzero_grid=round(len(znf525_binding_indices)/total * 50)
heatmap_g1 = fig.add_subplot(grid[0:nonzero_grid,20:70])
plots.plot_heatmap(phylop_447[znf525_binding_indices], cmap=custom_cmap.vlag_r_mpl,vmax=0.5,vmin=-0.5, matplot_axes = heatmap_g1)
heatmap_g1.xaxis.set_major_locator(MaxNLocator(integer=True))
heatmap_g1.set_yticks([])
heatmap_g1.set_title('phyloP_447 MER11 grouped by ZNF525 binding')
heatmap_g1.tick_params(axis='both', which='major', labelsize=8)
plots.plot_heatmap(znf808_forward[znf525_binding_indices], cmap = 'Blues',matplot_axes = heatmap_g1,vmin=0, vmax=0.1, opacity=0.5,)
plots.plot_heatmap(znf808_reverse[znf525_binding_indices], cmap = 'Reds',matplot_axes = heatmap_g1,vmin=0, vmax=0.1, opacity=0.5,)
plots.plot_heatmap(znf525_forward[znf525_binding_indices], cmap = 'Greens',matplot_axes = heatmap_g1,vmin=0, vmax=0.01, opacity=0.5,)
plots.plot_heatmap(znf525_reverse[znf525_binding_indices], cmap = 'Oranges',matplot_axes = heatmap_g1,vmin=0, vmax=0.01, opacity=0.5,)
############################
bar4.tick_params(axis='both', which='major', labelsize=8)
phylop_447_bind=phylop_447[znf525_binding_indices]
alignemnt_bind=alignment_filtered[znf525_binding_indices]
phylop_447_bind[alignemnt_bind == 0] = np.nan
phylop_447_nonbind=phylop_447[znf525_nonbinding_indices]
alignemnt_nonbind=alignment_filtered[znf525_nonbinding_indices]
phylop_447_nonbind[alignemnt_nonbind == 0] = np.nan
###################################
heatmap_g2 = fig.add_subplot(grid[nonzero_grid:50,20:70])
heatmap_g2.tick_params(axis='both', which='major', labelsize=8)
plots.plot_heatmap(phylop_447_nonbind, cmap=custom_cmap.vlag_r_mpl,vmax=0.5,vmin=-0.5, matplot_axes = heatmap_g2)
heatmap_g2.xaxis.set_major_locator(MaxNLocator(integer=True))
heatmap_g2.set_yticks([])
plots.plot_heatmap(znf808_forward[znf525_nonbinding_indices], cmap = 'Blues',matplot_axes = heatmap_g2,vmin=0, vmax=0.1, opacity=0.5,)
plots.plot_heatmap(znf808_reverse[znf525_nonbinding_indices], cmap = 'Reds',matplot_axes = heatmap_g2,vmin=0, vmax=0.1, opacity=0.5,)
plots.plot_heatmap(znf525_forward[znf525_nonbinding_indices], cmap = 'Greens',matplot_axes = heatmap_g2,vmin=0, vmax=0.01, opacity=0.5,)
plots.plot_heatmap(znf525_reverse[znf525_nonbinding_indices], cmap = 'Oranges',matplot_axes = heatmap_g2,vmin=0, vmax=0.01, opacity=0.5,)


bar = fig.add_subplot(grid[50:60,20:70])
bar.tick_params(axis='both', which='major', labelsize=8)
plots.plot_bar(phylop_447_bind, alignment=alignemnt_bind,mode='average', matplot_axes=bar,color='grey',opacity=0.5)
plots.plot_bar(phylop_447_nonbind, alignment=alignemnt_nonbind,mode='average', matplot_axes=bar,color ='grey',opacity=0.5)

bar2 = fig.add_subplot(grid[60:70,20:70])
bar2.tick_params(axis='both', which='major', labelsize=8)
plots.plot_bar(np.nanmedian(phylop_447_bind,axis=0)-np.nanmedian(phylop_447_nonbind,axis=0), alignment=alignment_filtered,mode='1D', matplot_axes=bar2)

bar3 = fig.add_subplot(grid[70:80,20:70])
bar3.tick_params(axis='both', which='major', labelsize=8)
plots.plot_bar(np.log(np.nanmedian(phylop_447_bind,axis=0)/np.nanmedian(phylop_447_nonbind,axis=0)), alignment=alignment_filtered,mode='1D', matplot_axes=bar3)

cbar = fig.add_subplot(grid[0:50, 70:72])
cbar.tick_params(axis='both', which='major', labelsize=8)
plots.plot_colorbar(cmap=custom_cmap.vlag_r_mpl,data=phylop_447, vmax=0.5,vmin=-0.5,matplot_axes=cbar, step = 0.1)
#plots.plot_colorbar(cmap='Purples',data=ap1,matplot_axes=cbar2, step = 0.1,vmin=0,vmax=1)
import matplotlib.ticker as ticker
def format_func(value, tick_number):
    if value < 0.1:
        return f"{value:.0e}"
    else:
        return f"{value:.2f}"
formatter = ticker.FuncFormatter(format_func)
heatmap_g1.margins(x=0, y=0)
heatmap_g2.margins(x=0, y=0)
######################################################age annotation

anno = fig.add_subplot(grid[0:50,18:20])
anno.tick_params(axis='both', which='major', labelsize=8)
anno.set_xticks([])
anno.set_yticks([])
cbar_anno =fig.add_subplot(grid[0:50,8:10])
cbar_anno.tick_params(axis='both', which='major', labelsize=8)
anno_sort = metadata_age.iloc[np.concatenate((znf525_binding_indices,znf525_nonbinding_indices))].te_age.fillna(0)
annot_plot = anno_sort.values
annot_uniq=anno_sort.sort_values().unique()

# Extract colors from the 'Blues' colormap
blues_cmap = plt.cm.Blues
num_colors = len(annot_uniq)
colors = [blues_cmap(i / (num_colors - 1)) for i in range(num_colors)]
# Create a colormap with a distinct color for each unique value
annot_cmap = mcolors.ListedColormap(colors)
# Add an extra boundary to center the ticks
boundaries = np.append(annot_uniq, annot_uniq[-1] + (annot_uniq[-1] - annot_uniq[-2]))
norm = mcolors.BoundaryNorm(boundaries, annot_cmap.N)

# Calculate tick positions to center them in the squares
tick_positions = (boundaries[:-1] + boundaries[1:]) / 2
# Create the plot
age_annot= anno.imshow(annot_plot.reshape(-1, 1),aspect = 'auto', cmap=annot_cmap, norm=norm)
#anno.set_title('age', fontsize='small')
# Create a colorbar with ticks at the calculated positions
anno_cbar=fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=annot_cmap),
             cax=cbar_anno, orientation='vertical', ticks= tick_positions)
anno_cbar.ax.set_yticklabels([str(val) for val in annot_uniq])
cbar_anno.set_title('(MYA)', fontsize='small')
cbar_anno.yaxis.set_ticks_position('left')
########################################################subgroup annotation
from matplotlib.colors import ListedColormap
anno2 = fig.add_subplot(grid[0:50,16:18])
anno2.tick_params(axis='both', which='major', labelsize=8)
anno2.set_xticks([])
cbar_anno2 =fig.add_subplot(grid[0:50,0:2])
cbar_anno2.tick_params(axis='both', which='major', labelsize=8)
metadata_age['subgroup'] = metadata_age['meta_id'].str.split('_').str[0]
subgroups = np.unique(metadata_age['subgroup'] )
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_annot=metadata_age['subgroup'].map(numerical_subgroup)[np.concatenate((znf525_binding_indices,znf525_nonbinding_indices))]

# Define a colormap (you can customize this)
colors = ['red', 'yellow', 'blue']  # Adjust these colors as needed
subgroup_cmap = ListedColormap(colors)
# Calculate boundaries for the categories
boundaries = np.arange(len(subgroups) + 1) - 0.5
s_norm = mcolors.BoundaryNorm(boundaries, subgroup_cmap.N)

age_annot= anno2.imshow(subgroup_annot.values.reshape(-1, 1),aspect = 'auto', cmap=subgroup_cmap)
#anno2.set_title('subgroup', fontsize='small')
anno_cbar2=fig.colorbar(mpl.cm.ScalarMappable(norm=s_norm, cmap=subgroup_cmap),
             cax=cbar_anno2, orientation='vertical', ticks= np.arange(len(subgroups)))
anno_cbar2.ax.set_yticklabels(['MER11A', 'MER11B', 'MER11C'])
cbar_anno2.set_title('subgroup', fontsize='small')
cbar_anno2.yaxis.set_ticks_position('left')
########################################################
import scipy
bar4 = fig.add_subplot(grid[80:90,20:70])

stat_v, p_value = scipy.stats.mannwhitneyu(phylop_447_bind,phylop_447_nonbind, axis =0,nan_policy='omit', alternative = 'less')
plots.plot_bar(-np.log(p_value), alignment=alignment_filtered,mode='1D', matplot_axes=bar4)
#########################################################
bar5 = fig.add_subplot(grid[90:95,20:70])
bar5.tick_params(axis='both', which='major', labelsize=8)
def min_max_scale(lst):
    if len(lst) == 1:
        return [0.0]  # or [0.5], depending on the desired output
    min_val = min(lst)
    max_val = max(lst)
    if min_val == max_val:
        return [0.0 for _ in lst]  # or [0.5 for _ in lst]
    return [(x - min_val) / (max_val - min_val) for x in lst]
plots.plot_bar(min_max_scale(znf808_minimum), alignment=alignment_filtered,mode='1D', matplot_axes=bar5,color='purple',opacity=0.5)
plots.plot_bar(min_max_scale(znf525_minimum), alignment=alignment_filtered,mode='1D', matplot_axes=bar5,color='brown',opacity=0.5)
cbar.set_title('phyloP', fontsize='small')
#cbar2.set_title('AP1', fontsize='small')
bar.set_title('phyloP', fontsize='small', x=1.06, y=0.4)
bar2.set_title('diff', fontsize='small', x=1.03, y=0.4)
bar3.set_title('log ratio', fontsize='small', x=1.07, y=0.4)
bar4.set_title('-log p-value\nranked test\nH1=less', fontsize='small', x=1.12, y=0.4)
bar5.set_title('KZFP binding', fontsize='small', x=1.1, y=0.4)
plt.show()
#%%