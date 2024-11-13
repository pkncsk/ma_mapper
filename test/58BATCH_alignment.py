#%% extract sequence
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import sequence_alignment
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
import os 
subfam_list= ['L1MB7','L1ME1', 'L1MEd',
       'L1MEf', 'L1MEg', 'L1PA16', 'L1PA3', 'L1PA4', 'L1PA5', 'L1PA7', 'L1PB1',
       'L2d', 'L3', 'MER20', 'MER3', 'MER58A', 'MER5A', 'MER5A1', 'MER5B',
       'MIR1_Amn', 'MLT1A0', 'MLT1B', 'MLT1C2', 'MLT1D', 'MLT1H', 'MLT1I',
       'MLT1J', 'MLT1K', 'MLT1L', 'MSTA', 'THE1B', 'THE1D', 'Tigger1']
#%%
for subfamily in subfam_list:
    print(f'process {subfamily}')
    fasta_file = f'{config.te_alignment_folder}/{subfamily}.fasta'
    
    if os.path.isfile(f'{fasta_file}.aligned'):
        print(f'{subfamily} already done')
        continue
    else:
        coord_file=sequence_alignment.extract_coord_from_repeatmasker_table(
            subfamily,
            repeatmasker_table = config.filtered_table, 
            #internal_id_table = f'{config.internal_id_folder}/{subfamily}.txt',
            save_to_file = True,
            output_filepath = f'{config.coord_internal_id_folder}/{subfamily}.txt')
        # alignment
        source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_fasta/hg38.fa'
        metadata = coord_file
        sequence_alignment.sequence_io(metadata,source_fasta,output_filepath =fasta_file, save_to_file= True,custom_id= False, custom_prefix = subfamily)
        if coord_file.shape[0]>1:
            sequence_alignment.mafft_align(fasta_file, nthread = 80, mafft_arg='--treeout --reorder --memsave ')

# %%
