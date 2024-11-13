#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_main as config
from ma_mapper import mapper
#%%
subfamily = 'MER11A'
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
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

# %%
coord_file = f'{config.coord_internal_id_folder}/{subfamily}.txt'
#%%
bam_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/KZFP-bam_hg38/znf808.sorted.bam'
#%%
bam_forward=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_forward')
bam_reverse=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_reverse')
bam_min=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_min')
#%%
from ma_mapper import plots
plots.all_overlay_plot(data = [bam_forward,bam_reverse], alignment= alignment_filtered,heatmap_annot=metadata_age.te_age, h_cmap=['Blues','Reds'], vlim = [[0,0.1],[0,0.1]], opacity = 0.5)
#%%