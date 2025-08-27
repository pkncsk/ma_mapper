#%% load package
from ma_mapper import mapper

alignment_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/alignment/MER11A.fasta.aligned'

alignment_matrix, alignment_coordinate  = mapper.parse_and_filter(alignment_file=alignment_filepath)
#%% MAF file
from ma_mapper import extract_maf
MAF_dir = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/'
output_matrix=mapper.map_and_overlay(alignment = alignment_filepath, data_file=MAF_dir,data_format='maf', separated_maf = True, count_arg='common_freq', target_species='hg38')

#%%
from ma_mapper import plots
#it is possible to change colormap
plots.plot(
    data=[output_matrix], 
    heatmap_color=["viridis"], 
    vlim =[[-0.5,0.5]],
    )
#%%
output_matrix
# %%
plots.plot(alignment=alignment_matrix,
    show_alignment=True,
    alignment_col = 'dna',
    heatmap_color=["viridis"], 
    vlim =[[-0.5,0.5]],
    )
# %%
