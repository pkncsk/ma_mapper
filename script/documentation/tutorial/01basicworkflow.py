#%% LOAD PACKAGE
from ma_mapper import mapper
from ma_mapper import plots
from ma_mapper import custom_cmap
#%% INITIAL PARAMETER
alignment_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/alignment/THE1C.fasta.aligned'
genomewide_data_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/annotation/homer_known_motif_hg38/AP-1(bZIP).bed'
#%% SIMPLE WORKFLOW
#extract genomewide data into data matrix using coordinates and alignment structure from the alignment file 
output_matrix=mapper.map_and_overlay(alignment_filepath, genomewide_data_filepath,data_format='bed')
#%% VISUALIZATION
#output matrix can be used for downstream analyses or visualization. this package also includes visualization helper function (matplotlib helper)
plots.plot(data = [output_matrix], heatmap_color=['Greens'], vlim = [[0,0.1]], opacity = 0.9)
#%% map_and_overlay() PIPELINE IN DETAIL
from ma_mapper import extract_bed
#extract TE coordinates from alignment fasta the header for the alignment should be in >NAME::chromosome:start-stop(strand) format
alignment_matrix, alignment_coordinate, filters  = mapper.parse_and_filter(alignment_file=alignment_filepath, preprocess_out=True)
#then use alignment coordinate to extract genomewide data
output_matrix=extract_bed.bed_io(coordinate_table=alignment_coordinate, bed=genomewide_data_filepath)
#filter coordinate table (exclude rows with content lower than a set threshold)
row_filter, col_filter = filters
alignment_coordinate_filtered=alignment_coordinate.iloc[row_filter]
#%%
#overlay/map extracted data matrix onto alignment matrix and filter rows and columns 
output_matrix_filtered = mapper.map_data(extracted_data=output_matrix, alignment_matrix= alignment_matrix, filter=filters)
#visualization
plots.plot(data = [output_matrix_filtered], heatmap_color=['Greens'], vlim = [[0,0.1]], opacity = 1)
#%% matrix extension
#it is possible to extend data extraction beyond the input coordinates by using extension_length argument. However the extended parts will be extracted as is, no alignment
phyloP_data_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/UCSC_phyloP_track/hg38.phyloP447way.bw'
phyloP_matrix=mapper.map_and_overlay(alignment_filepath, phyloP_data_filepath,data_format='bigwig')
plots.plot(data=[phyloP_matrix], heatmap_color=[custom_cmap.vlag_r_mpl], vlim =[[-0.5,0.5]])
phyloP_matrix_extended=mapper.map_and_overlay(alignment_filepath, phyloP_data_filepath,data_format='bigwig', extension_length=100)
plots.plot(data=[phyloP_matrix_extended], heatmap_color=[custom_cmap.vlag_r_mpl], vlim =[[-0.5,0.5]])
# %%
