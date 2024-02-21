#%%
from importlib import metadata
import sys
import pandas as pd
sys.path.append('/home/pc575/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import fetch_sequence
from ma_mapper import mafft_align
from ma_mapper import fetch_data
#%%

source_fasta = '/home/pc575/phd_project_development/data/hg38_fasta/hg38.fa'
metadata = '/home/pc575/phd_project_development/data/_ma_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11a_coord_with_id.txt'
fasta_file = '/home/pc575/phd_project_development/data/_ma_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/custom_id.fasta'
fetch_sequence.fetch_sequence(metadata,source_fasta,output_filepath =fasta_file, custom_id= True)

mafft_align.mafft_wrapper(fasta_file)

#bam_input= '/home/pc575/phd_project_development/data/functional_data/znf808.sorted.bam'
#metadata_input = '/home/pc575/phd_project_development/data/ma_mapper_output/mer11a_coord.txt'
#fetch_data.fetch_bam(metadata_input, bam_input)

# %%