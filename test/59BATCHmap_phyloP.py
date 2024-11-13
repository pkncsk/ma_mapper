#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import extract_bigwig, mapper, sequence_alignment
from ma_mapper import custom_cmap
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
import compress_pickle
import os 
#%%
subfam_list = config.subfam10k
#subfam_list = ['MER11A','MER11B','MER11C']
#%%
for subfamily in subfam_list:
    print(f'process {subfamily}')
    output_filepath = f'{config.phyloP_folder}/{subfamily}.lzma'
    #if os.path.isfile(output_filepath):
    #    print(f'{subfamily} already done')
    #    continue
    #else:
    alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
    if not os.path.isfile(alignment_file):
        print(f'alignment not found, try raw fasta')
        alignment_input =f'{config.te_alignment_folder}/{subfamily}.fasta'
        print(alignment_file)
    else:
        alignment_input = alignment_file
    alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_input)
    age_table = f'{config.te_age_folder}/{subfamily}.txt'
    age_df = pd.read_csv(age_table, sep='\t')
    internal_id_tbl = f'{config.internal_id_folder}/{subfamily}.txt'
    internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
    te_age_internal_id=internal_id_df.merge(age_df, on='internal_id')
    age_default_id = pd.DataFrame()
    age_default_id['internal_id'] = subfamily + '_' + te_age_internal_id.index.astype(str)
    age_default_id['te_age'] = te_age_internal_id['te_age']
    metadata_age = mapper.match_age_to_id_metadata(metadata_filtered, age_table=age_default_id)
    coord_file=mapper.extract_metadata_from_alignment(alignment_input)
    bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/zoonomia447/hg38.phyloP447way.bw'
    phylop_447=mapper.map_and_overlay(alignment_input, coord_file, bigwig_file, data_format='bigwig')
    compress_pickle.dump(phylop_447, output_filepath, compression="lzma")
# %%
#compress_pickle.load(output_filepath, compression="lzma")
# %%
