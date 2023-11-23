#%%
from importlib import metadata
import sys
import pandas as pd
sys.path.append('/home/pc575/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import fetch_sequence
#%%
metadata = '/home/pc575/phd_project_development/data/ma_mapper_output/mer11a_coord.txt'
source_fasta = '/home/pc575/phd_project_development/data/hg38_fasta/hg38.fa'
#%%
fetch_sequence.test_global()
#fetch_sequence.fetch_sequence(metadata,source_fasta)
# %%
