#%% LOAD PACKAGE
from ma_mapper import mapper
from ma_mapper import plots
#%% INITIAL PARAMETER
alignment_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/alignment/THE1C.fasta.aligned'
genomewide_data_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/annotation/homer_known_motif_hg38/AP-1(bZIP).bed'
#%% SIMPLE WORKFLOW
#extract genomewide data into data matrix using coordinates and alignment structure from the alignment file 
output_matrix=mapper.map_and_overlay(alignment_filepath, genomewide_data_filepath,data_format='bed', custom_id=True)
#%% VISUALIZATION
#output matrix can be used for downstream analyses or visualization. this package also includes visualization helper function (matplotlib helper)
plots.plot_experimental(data = [output_matrix], heatmap_color=['Greens'], vlim = [[0,0.1]], opacity = 1)
#%% map_and_overlay() PIPELINE IN DETAIL
from ma_mapper import extract_bed
#extract TE coordinates from alignment fasta the header for the alignment should be in >NAME::chromosome:start-stop(strand) format
alignment_matrix, alignment_coordinate, filters  = mapper.parse_and_filter(alignment_file=alignment_filepath, preprocess_out=True)
#then use alignment coordinate to extract genomewide data
output_matrix=extract_bed.bed_io(coordinate_table=alignment_coordinate, bed=genomewide_data_filepath, custom_id=True)
#filter coordinate table (exclude rows with content lower than a set threshold)
row_filter, col_filter = filters
alignment_coordinate_filtered=alignment_coordinate.iloc[row_filter]
#overlay/map extracted data matrix onto alignment matrix and filter rows and columns 
output_matrix_filtered = mapper.map_data(data_file=output_matrix, sorted_parsed_alignment= alignment_matrix, filter=filters)
#visualization
plots.plot_experimental(data = [output_matrix_filtered], heatmap_color=['Greens'], vlim = [[0,0.1]], opacity = 1)
#%%
