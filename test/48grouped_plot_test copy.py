from configparser import Interpolation
import sys
from tabnanny import verbose
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import mapper
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import plots
from ma_mapper import custom_cmap
import numpy as np
#%%
subfamily = ['THE1C']
alignment_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/alignment/'+subfamily[0]+'.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered, reference_table='/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/combined_age_div/combined_te_age_div.txt')
# %%
coord_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/coord_internal_id/'+subfamily[0]+'.txt'
#%%
#%%
bigwig_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_bigwig/241-mammalian-2020v2.bigWig'
phylop=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig')
#%%
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/AP-1(bZIP).bed'
ap1=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed')
#%% find highest
from scipy.signal import find_peaks
peaks, peaks_info = find_peaks(np.sum(ap1, axis = 0), width = 8)
highest_peak_index = peaks[np.argmax(np.sum(ap1, axis = 0)[peaks])]
#%%
nonzero_indices = np.where(ap1[:, highest_peak_index] != 0)[0]
zero_indices = np.where(ap1[:, highest_peak_index] == 0)[0]
#non_zero_indices = np.nonzero(np.any(ap1[:, highest_peak_indeces] != 0, axis=1))[0]
#zero_indices = np.nonzero(np.all(ap1[:, highest_peak_indeces] == 0, axis=1))[0]
#%%
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
grid = fig.add_gridspec(nrows = 100, ncols = 100, hspace=0)
total=len(zero_indices) + len(nonzero_indices)
nonzero_grid=round(len(nonzero_indices)/total * 50)

heatmap_g1 = fig.add_subplot(grid[0:nonzero_grid,20:70])
heatmap_g2 = fig.add_subplot(grid[nonzero_grid:50,20:70])
bar = fig.add_subplot(grid[50:60,20:70])
bar2 = fig.add_subplot(grid[60:70,20:70])
bar3 = fig.add_subplot(grid[70:80,20:70])
cbar = fig.add_subplot(grid[0:50, 70:72])
cbar2 = fig.add_subplot(grid[0:50, 77:79])
cbar_anno =fig.add_subplot(grid[0:50,10:12])
anno = fig.add_subplot(grid[0:50,18:20])
heatmap_g1.tick_params(axis='both', which='major', labelsize=8)
heatmap_g2.tick_params(axis='both', which='major', labelsize=8)
bar.tick_params(axis='both', which='major', labelsize=8)
cbar.tick_params(axis='both', which='major', labelsize=8)
cbar2.tick_params(axis='both', which='major', labelsize=8)
bar2.tick_params(axis='both', which='major', labelsize=8)
bar3.tick_params(axis='both', which='major', labelsize=8)
anno.tick_params(axis='both', which='major', labelsize=8)
cbar_anno.tick_params(axis='both', which='major', labelsize=8)
from ma_mapper import custom_cmap
plots.plot_heatmap(phylop[nonzero_indices], cmap=custom_cmap.vlag_r_mpl,vmax=0.5,vmin=-0.5, matplot_axes = heatmap_g1)
heatmap_g1.xaxis.set_major_locator(MaxNLocator(integer=True))
heatmap_g1.set_yticks([])
heatmap_g1.set_title('phyloP_241 THE1C grouped by AP1 motif')
plots.plot_heatmap(ap1[nonzero_indices], cmap = 'Greens',matplot_axes = heatmap_g1,vmin=0, vmax=1, opacity=0.5,)
plots.plot_heatmap(phylop[zero_indices], cmap=custom_cmap.vlag_r_mpl,vmax=0.5,vmin=-0.5, matplot_axes = heatmap_g2)
heatmap_g2.xaxis.set_major_locator(MaxNLocator(integer=True))
heatmap_g2.set_yticks([])
plots.plot_heatmap(ap1[zero_indices], cmap = 'Greens',matplot_axes = heatmap_g2,vmin=0, vmax=1, opacity=0.5,)
anno.set_xticks([])
plots.plot_bar(phylop[nonzero_indices], alignment=alignment_filtered[nonzero_indices],mode='average', matplot_axes=bar,color='blue',opacity=0.5)
plots.plot_bar(phylop[zero_indices], alignment=alignment_filtered[zero_indices],mode='average', matplot_axes=bar,color ='red',opacity=0.5)
plots.plot_bar(ap1, alignment=alignment_filtered,mode='average', matplot_axes=bar,color='green',opacity=0.5, ylim=bar.get_ylim())

plots.plot_bar(np.nanmean(phylop[nonzero_indices],axis=0)-np.nanmean(phylop[zero_indices],axis=0), alignment=alignment_filtered,mode='1D', matplot_axes=bar2)
plots.plot_bar(ap1, alignment=alignment_filtered,mode='average', matplot_axes=bar2,color='green',opacity=0.2, ylim=bar2.get_ylim())

plots.plot_bar(np.log(np.nanmean(phylop[nonzero_indices],axis=0)/np.nanmean(phylop[zero_indices],axis=0)), alignment=alignment_filtered,mode='1D', matplot_axes=bar3)
plots.plot_bar(ap1, alignment=alignment_filtered,mode='average', matplot_axes=bar3,color='green',opacity=0.2, ylim=bar3.get_ylim())

plots.plot_colorbar(cmap=custom_cmap.vlag_r_mpl,data=phylop, vmax=0.5,vmin=-0.5,matplot_axes=cbar, step = 0.1)
plots.plot_colorbar(cmap='Greens',data=ap1,matplot_axes=cbar2, step = 0.1,vmin=0,vmax=1)
import matplotlib.ticker as ticker
def format_func(value, tick_number):
    if value < 0.1:
        return f"{value:.0e}"
    else:
        return f"{value:.2f}"
formatter = ticker.FuncFormatter(format_func)

heatmap_g1.margins(x=0, y=0)
heatmap_g2.margins(x=0, y=0)
anno_sort = metadata_age.iloc[np.concatenate((nonzero_indices,zero_indices))].te_age.fillna(0)
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
norm = mcolors.BoundaryNorm(boundaries, annot_cmap.N)
# Calculate tick positions to center them in the squares
tick_positions = (boundaries[:-1] + boundaries[1:]) / 2
# Create the plot
age_annot= anno.imshow(annot_plot.reshape(-1, 1),aspect = 'auto', cmap=annot_cmap, norm=norm)
anno.set_title('age', fontsize='small')
# Create a colorbar with ticks at the calculated positions
anno_cbar=fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=annot_cmap),
             cax=cbar_anno, orientation='vertical', ticks= tick_positions)
anno_cbar.ax.set_yticklabels([str(val) for val in annot_uniq])
cbar_anno.set_title('(MYA)', fontsize='small')
cbar_anno.yaxis.set_ticks_position('left')
import scipy
bar4 = fig.add_subplot(grid[80:90,20:70])
bar4.tick_params(axis='both', which='major', labelsize=8)
stat_v, p_value = scipy.stats.mannwhitneyu(phylop[nonzero_indices],phylop[zero_indices], axis =0,nan_policy='omit',alternative='greater')
plots.plot_bar(-np.log(p_value), alignment=alignment_filtered,mode='1D', matplot_axes=bar4)
plots.plot_bar(ap1, alignment=alignment_filtered,mode='average', matplot_axes=bar4,color='green',opacity=0.2, ylim=bar4.get_ylim())
cbar.set_title('phyloP', fontsize='small')
cbar2.set_title('AP1', fontsize='small')
bar.set_title('phyloP', fontsize='small', x=1.06, y=0.4)
bar2.set_title('diff', fontsize='small', x=1.03, y=0.4)
bar3.set_title('log ratio', fontsize='small', x=1.07, y=0.4)
bar4.set_title('-log p-value\nranked test\nH1=greater than', fontsize='small', x=1.12, y=0.4)
plt.show()

# %%
