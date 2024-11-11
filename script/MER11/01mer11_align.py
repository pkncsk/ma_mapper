#%% extract sequence
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import sequence_alignment
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_main as config
subfamily = ['MER11A','MER11B','MER11C']
coord_file=sequence_alignment.extract_coord_from_repeatmasker_table(
    subfamily,
    repeatmasker_table = config.filtered_table, 
    internal_id_table = [
    '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker/internal_id/MER11A.txt',
    '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker/internal_id/MER11B.txt',
    '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker/internal_id/MER11C.txt',],
    save_to_file = False,
    output_filepath = config.coord_internal_id_folder+'/MER11.txt')
#%% alignment
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from Bio.SeqRecord import SeqRecord
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import sequence_alignment

subfamily = 'MER11'
source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_fasta/hg38.fa'
metadata = coord_file
fasta_file = f'{config.te_alignment_folder}/{subfamily}.fasta'
#%%
sequence_alignment.sequence_io(metadata,source_fasta,output_filepath =fasta_file, save_to_file= True,custom_id= False, custom_prefix=subfamily)
#%%
sequence_alignment.mafft_align(fasta_file, nthread = 6)
#%% parse alignment
import sys
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_main as config
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
subfamily = 'MER11'
from ma_mapper import mapper
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
#%%
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file, col_threshold = 0, col_content_threshold = 0, row_threshold = 0)
age_table_list = [f'{config.te_age_folder}/MER11A.txt',
             f'{config.te_age_folder}/MER11B.txt',
             f'{config.te_age_folder}/MER11C.txt']
age_df_list = []
for age_tbl in age_table_list:
    age_df_list.append(pd.read_csv(age_tbl, sep='\t'))
age_df=pd.concat(age_df_list)
internal_id_tbl = f'{config.internal_id_folder}/{subfamily}.txt'
internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
te_age_internal_id=internal_id_df.merge(age_df, on='internal_id')
age_default_id = pd.DataFrame()
age_default_id['internal_id'] = subfamily + '_' + te_age_internal_id.index.astype(str)
age_default_id['te_age'] = te_age_internal_id['te_age']
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered, age_table= age_default_id)
#%%
from ma_mapper import plots
plots.plot_experimental(alignment=alignment_filtered, alignment_col='dna', show_alignment_colbar=True,interpolation ='nearest')
#%%
alignment_extend, metadata_filtered= mapper.parse_and_filter(alignment_file, extension_length=500)
#%%
plots.all_overlay_plot(alignment=alignment_extend, alignment_col='dna', interpolation ='nearest')
# %%
