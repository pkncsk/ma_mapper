#%% LOAD PACKAGE
from ma_mapper import mapper
from ma_mapper import plots
#%% INITIAL PARAMETER
alignment_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/alignment/THE1C.fasta.aligned'
genomewide_data_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/annotation/homer_known_motif_hg38/AP-1(bZIP).bed'
#extract genomewide data into data matrix using coordinates and alignment structure from the alignment file 
output_matrix=mapper.map_and_overlay(alignment_filepath, genomewide_data_filepath,data_format='bed', custom_id=True)
#%% DOWNSTREAM ANALYSES
# as stated in the introduction, ma_mapper output is a numerical matrix of genome-wide data that has the same gap-structure like sequence alignment
# to properly analyze the extracted data, ma_mapper offer a function that would consider the gaps structure in each column using alignment matrix as reference 
alignment_matrix, coordinate_table=mapper.parse_and_filter(alignment_filepath)
mean_output=mapper.normalise(alignment=alignment_matrix, mapped_data=output_matrix, method = 'average')
coverage_output=mapper.normalise(alignment=alignment_matrix, mapped_data=output_matrix, method = 'perc_coverage')
median_output=mapper.normalise(alignment=alignment_matrix, mapped_data=output_matrix, method = 'median')
#%% VISUALIZATION
# the plot function can also helps plotting a 2D plot
plots.plot_experimental(aggregated_data=[mean_output], aggregated=True,heatmap=False,agg_colset=['grey'], agg_ylabel=['average value'],agg_ylabel_fs=8, agg_xlabel='position (bp)', agg_xlabel_fs=8)
#%% plots can be stacked 
plots.plot_experimental(aggregated_data=[mean_output, median_output, coverage_output], aggregated=True,heatmap=False,agg_colset=['grey','grey','red'], agg_ylabel=['average','median','coverage'],agg_ylabel_fs=8, agg_xlabel='position (bp)', agg_xlabel_fs=8)
#%% CUSTOM ANALYSIS
# since the plot function accepts 1d array as long as the length matches, a user can make their own analysis and plot the result using plot
# for example