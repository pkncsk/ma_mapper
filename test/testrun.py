#%%
from importlib import metadata
import sys
import pandas as pd
sys.path.append('/home/pc575/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import fetch_sequence
from ma_mapper import mafft_align
#%%
def main():
    #source_fasta = '/home/pc575/phd_project_development/data/hg38_fasta/hg38.fa'
    #metadata = '/home/pc575/phd_project_development/data/ma_mapper_output/mer11a_coord.txt'
    #fetch_sequence.fetch_sequence(metadata,source_fasta, nthread= 2)
    mer11_fasta = '/home/pc575/phd_project_development/data/ma_mapper_output/seqrecords.fasta'
    mafft_align.mafft_wrapper(mer11_fasta)
#%%
if __name__ == '__main__':
    main()
# %%
