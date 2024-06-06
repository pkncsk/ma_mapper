#%% extract sequence
from math import e
import re
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
subfamily = ['MER11A','MER11B','MER11C']
from ma_mapper import sequence_alignment
#coord_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/coord_internal_id/'+subfamily[0]+'.txt'
coord_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11/MER11A.txt'
sequence_alignment.extract_subfamily_coord(subfamily, output_filepath= coord_file, species_reference='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/combined_age_div/combined_age_and_div.txt')
#%%
import sys
import os
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
internal_id_table = None
repeatmasker_table = config.filtered_table
subfamilies = ['MER11A','MER11B','MER11C']
for idx,subfamily in enumerate(subfamilies):
    if isinstance(repeatmasker_table, str):
        repeatmasker_table=pd.read_csv(repeatmasker_table, sep='\t', index_col = 0, header = 0)
    else:
        repeatmasker_table=repeatmasker_table
    if internal_id_table is not None:
        if isinstance(internal_id_table, list):
            internal_id_tbl = internal_id_table[idx]
        else:
            internal_id_tbl = internal_id_table
        if isinstance(internal_id_tbl, str):
            if os.path.isfile(internal_id_tbl):
                internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
        elif isinstance(internal_id_tbl, pd.DataFrame):
            internal_id_df = internal_id_tbl
        subfam_w_internal_id=pd.merge(repeatmasker_table, internal_id_df, left_index = True, right_on = 'rmsk_index')
        subfam_coords = subfam_w_internal_id[['genoName','genoStart','genoEnd','internal_id']]
        subfam_coords['score'] = 10
        subfam_coords['strand'] = subfam_w_internal_id.strand
    else:
        subfam_table=repeatmasker_table[repeatmasker_table.repName == subfamilies]
#%% alignment
import sys
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import sequence_alignment
subfamily = ['THE1C']
source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_fasta/hg38.fa'
metadata = coord_file
fasta_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/alignment/'+subfamily[0]+'.fasta'
#%%
sequence_alignment.fetch_sequence(metadata,source_fasta,output_filepath =fasta_file, save_to_file= True,custom_id= True)
sequence_alignment.mafft_align(fasta_file, nthread = 40)
#%% parse alignment
import sys
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
subfamily = ['THE1C']
from ma_mapper import mapper
alignment_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/alignment/'+subfamily[0]+'.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered)
#%%
from ma_mapper import plot
plot.overlay_plot([alignment_filtered], metadata= metadata_age)
#%%
alignment_extend, metadata_filtered= mapper.parse_and_filter(alignment_file, extension_length=500)
plot.overlay_plot([alignment_extend], metadata= metadata_age)
#%%
