#%% LOAD PACKAGE
from ma_mapper import mapper
from ma_mapper import plots
from ma_mapper import custom_cmap
#%% PLOT CONFIGURATION
#the plots module is a matplotlib wrapper that offer a certain degree of quick configuration for mapped alignment visualization

#plot_experimental() is the wrapper of functions in this module. At basic level, this function should be enough for data plotting and customization
#%% plot preparation (not shown)
alignment_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/alignment/THE1C.fasta.aligned'
AP1_motif_data_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/annotation/homer_known_motif_hg38/AP-1(bZIP).bed'
NFkB_motif_data_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/annotation/homer_known_motif_hg38/NFkB-p65(RHD).bed'
phyloP_data_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/UCSC_phyloP_track/hg38.phyloP447way.bw'

#extract genomewide data into data matrix using coordinates and alignment structure from the alignment file 
ap1_matrix=mapper.map_and_overlay(alignment_filepath, AP1_motif_data_filepath,data_format='bed', custom_id=True)
nfkb_matrix=mapper.map_and_overlay(alignment_filepath, NFkB_motif_data_filepath,data_format='bed', custom_id=True)
phyloP_matrix = mapper.map_and_overlay(alignment_filepath, phyloP_data_filepath, data_format='bigwig', custom_id=True)
#%%heatmap customization
#as the output from mapper module is matrix, one of the best way to visualize it is heatmap, so this module would play around heatmap plot
#here is the most basic form of heatmap
plots.plot_experimental(
    data=[phyloP_matrix], #data matrix
    heatmap_color=[custom_cmap.vlag_r_mpl], #colormap for the heatmap
    vlim =[[-0.5,0.5]]) #data value caps
#%%
#it is possible to change colormap
plots.plot_experimental(
    data=[phyloP_matrix], 
    heatmap_color=["viridis"], 
    vlim =[[-0.5,0.5]],
    )
#%% there might be some data that need to be shown with alignment gap for visual/context clarity, there are specific functions to handle alignment plot
alignment_matrix, coordinate_table=mapper.parse_and_filter(alignment_filepath)
plots.plot_experimental(
    show_alignment=True,
    alignment=alignment_matrix, 
    alignment_col='dna', # color coding for dna
    )
#%%
# in some context, we only need to see the gaps
plots.plot_experimental(
    show_alignment=True,
    alignment=alignment_matrix, 
    alignment_col='nulc_white', 
    )
#%% we can also add colorbars for both heatmap and alignment
plots.plot_experimental(
    data=[phyloP_matrix], 
    heatmap_color=["viridis"], 
    vlim =[[-0.5,0.5]],
    colorbar = True, #enable colorbar
    colorbar_steps = [0.1], #control colorbar scale
    )
#%%
plots.plot_experimental(
    show_alignment=True,
    alignment=alignment_matrix, 
    alignment_col='dna',
    colorbar = True,
    show_alignment_colbar=True 
    )
#%% similar to data overlay, we can overlay data over alignment
plots.plot_experimental(
    data=[ap1_matrix], 
    heatmap_color=["Blues"], 
    vlim =[[0,5]],
    opacity=0.99, #transparency value of heatmap layer
    transparency_mode = 'gradient', 
    show_alignment=True,
    alignment=alignment_matrix, 
    alignment_col='nulc_white', 
    )
#%%
#title, x-, y- axis label can be added 
plots.plot_experimental(
    data=[phyloP_matrix], 
    heatmap_color=[custom_cmap.vlag_r_mpl], 
    vlim =[[-0.5,0.5]],
    heatmap_title=['example plot'],
    heatmap_title_fs=6, #font size
    heatmap_ylabel='seqeunces',
    heatmap_ylabel_fs=14,
    heatmap_xlabel='position (bp)',
    heatmap_xlabel_fs=16,
    ) 
# %%
# also a certain area of the plot an be highlighted
plots.plot_experimental(
    data=[phyloP_matrix], 
    heatmap_color=[custom_cmap.vlag_r_mpl], 
    vlim =[[-0.5,0.5]],
    heatmap_yhighlight= [[100,200]],#cooridnates
    heatmap_yhighlight_col= ['red'], #highlight color
    heatmap_yhighlight_alpha = [0.2], #transparent value
    heatmap_xhighlight = [[100,200],[6000,800]],
    heatmap_xhighlight_col = ['blue','green'],
    heatmap_xhighlight_alpha = [0.5,0.5],
    )
#%% aggregated plot
# this part was intended to show aggregated data from the main plot by normalizing the matrix by column and plot a resulting array as a bar chart
##plot preparation (not shown)
mean_phyloP=mapper.normalise(alignment=alignment_matrix, mapped_data=phyloP_matrix, method = 'average')
#%%
#here is the most basic aggregated plot
plots.plot_experimental(
    heatmap=False,
    aggregated_data=[mean_phyloP], 
    aggregated=True,
    agg_colset=['grey'])
# %% same as the heatmap, it is possible to change plot color
plots.plot_experimental(
    heatmap=False,
    aggregated_data=[mean_phyloP], 
    aggregated=True,
    agg_colset=['red'])
#%%
#title, x-, y- axis label can be added 
plots.plot_experimental(
    heatmap=False,
    aggregated_data=[mean_phyloP], 
    aggregated=True,
    agg_colset=['grey'],
    agg_titles=['example 1d plot'],
    agg_titles_fs=12,
    agg_xlabel='position (bp)',
    agg_xlabel_fs=4,
    agg_ylabel=['phyloP'],
    agg_ylabel_fs=8)

# %%
# also a certain area of the plot an be highlighted
plots.plot_experimental(
    heatmap=False,
    aggregated_data=[mean_phyloP], 
    aggregated=True,
    agg_colset=['grey'],
    agg_xhighlight=[[-0.5,0]],
    agg_xhighlight_col=['red'],
    agg_xhighlight_alpha=[0.5],
    agg_yhighlight=[[50,100],[200,210]],
    agg_yhighlight_col=['blue','green'],
    agg_yhighlight_alpha=[0.1,0.9])
#%% annotation configuration
#the annotation and annotatin colorbar can be modified to a certain degree
##plot preparation (not shown)
#genomewide data can be used to sort data matrix
#calculate AP-1 coverage
ap1_coverage_array=mapper.normalise(alignment=alignment_matrix, mapped_data=ap1_matrix, method='perc_coverage')
#since the coordinate table was filtered in parse and filter function, reset index to match positions in output matrices
coordinate_table = coordinate_table.reset_index()
#find peaks of coverage
import scipy
peaks, _ = scipy.signal.find_peaks(ap1_coverage_array, width = 6)
import numpy as np
index_of_rows_with_ap1 = np.unique(np.where(ap1_matrix[:, peaks] != 0)[0])
index_of_rows_without_ap1=list(set(np.arange(ap1_matrix.shape[0])) - set(index_of_rows_with_ap1))
coordinate_table['AP1_motif'] = 1
coordinate_table.loc[coordinate_table.index.isin(index_of_rows_with_ap1), 'AP1_motif'] = 0
#extract annotation from metadata
ap1_motif_annotation=coordinate_table['AP1_motif']
# %% basic annotation configuration
plots.plot_experimental(
    data=[phyloP_matrix], 
    heatmap_color=[custom_cmap.vlag_r_mpl], 
    vlim =[[-0.5,0.5]], 
    opacity=0.99, 
    transparency_mode = 'gradient', 
    annotation = True,  #enable annotation
    anno_col=[['blue','white']], #annotation colorset
    annotation_data= [ap1_motif_annotation], #annotation array
    anno_cbar_label=[['TE with AP1 motif', 'TE without AP1 motif']]) #annotation label
# %% we can adjust label and title
plots.plot_experimental(
    data=[phyloP_matrix], 
    heatmap_color=[custom_cmap.vlag_r_mpl], 
    vlim =[[-0.5,0.5]], 
    opacity=0.99, 
    transparency_mode = 'gradient', 
    annotation = True,  #enable annotation
    anno_col=[['blue','white']], #annotation colorset
    annotation_data= [ap1_motif_annotation], #annotation array
    anno_cbar_label=[['TE with AP1 motif', 'TE without AP1 motif']], #annotation label
    anno_ylabel='example',
    anno_ylabel_fs=12,
    anno_cbar = True,
    anno_cbar_title=['example cbar'],
)
# %% 
