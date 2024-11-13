#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
from ma_mapper import mapper
from ma_mapper import custom_cmap
#%%
subfamily = 'MER11'
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file, col_threshold = 0.5, row_threshold = 0.5)
metadata_filtered['len'] = metadata_filtered.end.astype(int) - metadata_filtered.start.astype(int)
print(alignment_filtered.shape)

#%%
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
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/zoonomia447/hg38.phyloP447way.bw'
phylop_447=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', pf_col_threshold = 0.5, pf_row_threshold = 0.5)
# %%
from ma_mapper import plots
plots.all_overlay_plot(data = [phylop_447], alignment= alignment_filtered,heatmap_annot=metadata_age.te_age, h_cmap=[custom_cmap.vlag_r_mpl], vlim = [[-0.5,0.5]], opacity = 0.5,)
#%%