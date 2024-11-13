#%% extract sequence
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')

from ma_mapper import sequence_alignment
subfamily = 'protein_coding_sequences_sliced'
coord_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_main/gene_ref/protein_coding_sequences_sliced.txt'
#%% alignment
import sys
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import sequence_alignment
subfamily = 'protein_coding_sequences_sliced'
source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_fasta/hg38.fa'
metadata = coord_file
fasta_file = f'/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_main/alignment/{subfamily}.fasta'
#%%
sequence_alignment.sequence_io(metadata,source_fasta,output_filepath =fasta_file, save_to_file= True,custom_id= True)
#%%
#sequence_alignment.mafft_align(fasta_file, nthread = 40)
#%% parse alignment
import sys
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
subfamily = 'protein_coding_sequences_sliced'
from ma_mapper import mapper
alignment_file = f'/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_main/alignment/{subfamily}.fasta'
alignment_extended, metadata_filtered= mapper.parse_and_filter(alignment_file,custom_id=True, extension_length=500, source_fasta = source_fasta)
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(heatmap=True,
    show_alignment=True,
    alignment=alignment_extended, 
    heatmap_mode='overlay', 
    alignment_col='dna',
    show_alignment_colbar=True,
    colorbar=True,
    opacity = 0.5, 
    agg_major_tick=200,
    #figsize=[100,30],

    )
#%%