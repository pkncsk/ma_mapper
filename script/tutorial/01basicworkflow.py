#%% LOAD PACKAGE
from ma_mapper import mapper
from ma_mapper import plots
#%% INITIAL PARAMETER
alignment_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/alignment/THE1C.fasta.aligned'
genomewide_data_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/annotation/homer_known_motif_hg38/AP-1(bZIP).bed'
#%% BASIC WORKFLOW
#extract coordinates from alignment file (can be created by other means, only need to be a BED file)
coordinate_table= mapper.extract_coordinate_from_alignment(alignment_filepath)
#extract genomewide data into data matrix using 1) coordinate table for data extraction, 2) alignment file for container array reference
output_matrix=mapper.map_and_overlay(alignment_filepath, genomewide_data_filepath,data_format='bed', coordinate_file = coordinate_table, custom_id=True)
#%% VISUALIZATION
#output matrix can be used for downstream analyses or visualization
plots.plot_experimental(data = [output_matrix], heatmap_color=['Greens'], vlim = [[0,0.1]], opacity = 1)
#%% MORE DETAILED STEP-BY-STEP (mapper.map_and_overlay() pipeline)
#prepare TE coordinate file (BED format, can be extracted from alignment or use obtained by other means)
coordinate_table= mapper.extract_metadata_from_alignment(alignment_filepath)
#parse the alignment file into a numerical matrix and create filter based on alignment content, coordinate extrated from alignment will be used as order reference 
alignment_numerical_matrix, metadata_table= mapper.parse_and_filter(alignment_filepath)


# %%
alignment_matrix, alignment_coordinate, filters  = mapper.parse_and_filter(alignment_file=alignment_filepath, preprocess_out=True)
# %%
from ma_mapper import extract_bed
output_matrix=extract_bed.bed_io(coordinate_table=alignment_coordinate, bed=genomewide_data_filepath, custom_id=True)
# %%
import numpy as np
source_order = alignment_coordinate.iloc[:,3].unique()
sorted_order = []
for _, row in alignment_coordinate.iterrows():
    sorted_order.append(np.where(source_order == row['name'])[0][0])
# %%
