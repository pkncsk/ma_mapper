#%%
from importlib import metadata
import sys
import pandas as pd
sys.path.append('/home/pc575/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import fetch_sequence
from ma_mapper import mafft_align
from ma_mapper import mapper
#%%
def main():
    #source_fasta = '/home/pc575/phd_project_development/data/hg38_fasta/hg38.fa'
    #metadata = '/home/pc575/phd_project_development/data/ma_mapper_output/mer11a_coord.txt'
    #fetch_sequence.fetch_sequence(metadata,source_fasta, nthread= 2)
    #mer11_fasta = '/home/pc575/phd_project_development/data/ma_mapper_output/seqrecords.fasta'
    #mafft_align.mafft_wrapper(mer11_fasta)
        
    bam_input= '/home/pc575/phd_project_development/data/functional_data/znf808.sorted.bam'
    metadata_input = '/home/pc575/phd_project_development/data/ma_mapper_output/mer11a_coord.txt'
    mapper.fetch_bam(metadata_input, bam_input)
#%%
if __name__ == '__main__':
    main()
# %%