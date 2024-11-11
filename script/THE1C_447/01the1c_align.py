#%% extract sequence
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import sequence_alignment
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
subfamily = 'THE1C'
coord_file=sequence_alignment.extract_coord_from_repeatmasker_table(
    subfamily,
    repeatmasker_table = config.filtered_table, 
    #internal_id_table = f'{config.internal_id_folder}/{subfamily}.txt',
    save_to_file = True,
    output_filepath = f'{config.coord_internal_id_folder}/{subfamily}.txt')
#%% alignment
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import sequence_alignment

subfamily = 'THE1C'
source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_fasta/hg38.fa'
metadata = coord_file
fasta_file = f'{config.te_alignment_folder}/{subfamily}.fasta'
#%%
sequence_alignment.sequence_io(metadata,source_fasta,output_filepath =fasta_file, save_to_file= True,custom_id= False, custom_prefix = subfamily)
#%%
sequence_alignment.mafft_align(fasta_file, nthread = 6)
#%% parse alignment
import sys
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
subfamily = 'THE1C'
from ma_mapper import mapper
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
import pandas as pd
#%%
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
age_table = f'{config.te_age_folder}/{subfamily}.txt'
age_df = pd.read_csv(age_table, sep='\t')
internal_id_tbl = f'{config.internal_id_folder}/{subfamily}.txt'
internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
te_age_internal_id=internal_id_df.merge(age_df, on='internal_id')
age_default_id = pd.DataFrame()
age_default_id['internal_id'] = subfamily + '_' + te_age_internal_id.index.astype(str)
age_default_id['te_age'] = te_age_internal_id['te_age']
#%%
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered, age_table= age_default_id)
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
plots.plot_experimental(alignment=alignment_filtered, alignment_col='dna',show_alignment=True, show_alignment_colbar=True, colorbar=True)
#%%
alignment_extend, metadata_filtered= mapper.parse_and_filter(alignment_file, extension_length=500, source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_fasta/hg38.fa')
#%%
plots.plot_experimental(alignment=alignment_extend, alignment_col='dna', show_alignment=True, show_alignment_colbar=True, colorbar=True,figsize=[100,20])
# %%
