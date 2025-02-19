#%% LOAD PACKAGE
from ma_mapper import mapper
from ma_mapper import plots
from ma_mapper import custom_cmap
import pandas as pd
#%% INITIAL PARAMETER
alignment_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/alignment/THE1C.fasta.aligned'
TF_motif_data_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/annotation/homer_known_motif_hg38/AP-1(bZIP).bed'
phyloP_data_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/UCSC_phyloP_track/hg38.phyloP447way.bw'
#%%
#extract genomewide data into data matrix using coordinates and alignment structure from the alignment file 
ap1_matrix=mapper.map_and_overlay(alignment_filepath, TF_motif_data_filepath,data_format='bed', custom_id=True)
#%%
phyloP_matrix = mapper.map_and_overlay(alignment_filepath, phyloP_data_filepath, data_format='bigwig', custom_id=True)
#%% DATA OVERLAY
# Since different kinds of data are overlaid on the same TE alignment, it is possible to overlay their plot to gain visual insight about the data 
plots.plot_experimental(data=[phyloP_matrix,ap1_matrix], heatmap_color=[custom_cmap.vlag_r_mpl,"Greens"], vlim =[[-0.5,0.5],[0,10]], opacity=0.5)
#%% PLOT ANNOTATION
#data matrix extracted from genome-wide data of TE can also be annotated with different data that can link to TEs which can be included to the data matrix using coordinate table as a metadata table.