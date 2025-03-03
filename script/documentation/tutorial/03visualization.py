#%% LOAD PACKAGE
from ma_mapper import mapper
from ma_mapper import plots
from ma_mapper import custom_cmap
import pandas as pd
#%% INITIAL PARAMETER
alignment_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/alignment/THE1C.fasta.aligned'
AP1_motif_data_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/annotation/homer_known_motif_hg38/AP-1(bZIP).bed'
NFkB_motif_data_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/annotation/homer_known_motif_hg38/NFkB-p65(RHD).bed'
phyloP_data_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/UCSC_phyloP_track/hg38.phyloP447way.bw'
#%%
#extract genomewide data into data matrix using coordinates and alignment structure from the alignment file 
ap1_matrix=mapper.map_and_overlay(alignment_filepath, AP1_motif_data_filepath,data_format='bed')
nfkb_matrix=mapper.map_and_overlay(alignment_filepath, NFkB_motif_data_filepath,data_format='bed')
#%% DATA OVERLAY
# Since different kinds of data are overlaid on the same TE alignment, it is possible to overlay their plot to gain visual insight about the data 
plots.plot(data=[nfkb_matrix,ap1_matrix], heatmap_color=["Blues","Greens"], vlim =[[0,10],[0,10]], opacity=0.99)
#%%  PLOT SORTING
#genomewide data can be used to sort data matrix
#calculate AP-1 coverage
alignment_matrix, coordinate_table=mapper.parse_and_filter(alignment_filepath)
ap1_coverage_array=mapper.normalise(alignment_matrix=alignment_matrix, data_matrix=ap1_matrix, method='perc_coverage')
#since the coordinate table was filtered in parse and filter function, reset index to match positions in output matrices
coordinate_table = coordinate_table.reset_index()
#find peaks of coverage
import scipy
peaks, _ = scipy.signal.find_peaks(ap1_coverage_array, width = 6)
import numpy as np
index_of_rows_with_ap1 = np.unique(np.where(ap1_matrix[:, peaks] != 0)[0])
index_of_rows_without_ap1=list(set(np.arange(ap1_matrix.shape[0])) - set(index_of_rows_with_ap1))
index_sorted=np.concatenate((index_of_rows_with_ap1, index_of_rows_without_ap1))
#the sorted index can be used to rearrange matrix of mapped genomewide data
ap1_matrix_sorted = ap1_matrix[index_sorted]
nfkb_matrix_sorted = nfkb_matrix[index_sorted]
#%%
plots.plot(data=[nfkb_matrix_sorted,ap1_matrix_sorted], heatmap_color=["Blues","Greens"], vlim =[[0,10],[0,10]], opacity=0.99)
#%% PLOT ANNOTATION
#some data such as phyloP can use visual aid from TF motif annotation
phyloP_matrix = mapper.map_and_overlay(alignment_filepath, phyloP_data_filepath, data_format='bigwig')
#data matrix extracted from genome-wide data of TE can also be annotated using coordinate table as a metadata table.
coordinate_table['AP1_motif'] = 1
coordinate_table.loc[coordinate_table.index.isin(index_of_rows_with_ap1), 'AP1_motif'] = 0
#%%
#extract annotation from metadata
ap1_motif_annotation=coordinate_table['AP1_motif']
plots.plot(data=[phyloP_matrix], heatmap_color=[custom_cmap.vlag_r_mpl], vlim =[[-0.5,0.5]], opacity=0.99, annotation = True, anno_col=[['blue','white']], annotation_data= [ap1_motif_annotation], anno_cbar_label=[['TE with AP1 motif', 'TE without AP1 motif']])
#%% ANNOTATION SORTING
phyloP_matrix_sorted = phyloP_matrix[index_sorted]
plots.plot(data=[phyloP_matrix], heatmap_color=[custom_cmap.vlag_r_mpl], vlim =[[-0.5,0.5]], opacity=0.99)
plots.plot(data=[phyloP_matrix_sorted], heatmap_color=[custom_cmap.vlag_r_mpl], vlim =[[-0.5,0.5]], opacity=0.99)
#sort rows in metadata using row order from coverage sort earlier
coordinate_table_sorted=coordinate_table.iloc[index_sorted]
#now the metadata table has the same order as the sorted alignment, extract annotation
ap1_motif_annotation_sorted=coordinate_table_sorted['AP1_motif']
plots.plot(data=[phyloP_matrix_sorted], heatmap_color=[custom_cmap.vlag_r_mpl], vlim =[[-0.5,0.5]], opacity=0.99, annotation = True, anno_col=[['blue','white']], annotation_data= [ap1_motif_annotation_sorted], anno_cbar_label=[['TE with AP1 motif', 'TE without AP1 motif']])

# %%
