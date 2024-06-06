#%%
import pandas as pd
import sys
import os
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/stable')
import config_hg38_repeatlib as config
from concurrent.futures import ProcessPoolExecutor
#%%
def make_te_age_bed(subfamily):

#    subfamily = 'THE1C'
    subfamily_filename = subfamily.replace('/','_')
    print(f'process: {subfamily}')
    #subfamily = 'MER5C'
    #coord_internal_id_folder = config.coord_internal_id_folder
    #te_coord_tbl = f'{coord_internal_id_folder}/{subfamily_filename}.txt'
    #te_coord_df = pd.read_csv(te_coord_tbl, sep='\t')
    filtered_table = config.filtered_table
    internal_id_folder = config.internal_id_folder
    internal_id_tbl = f'{internal_id_folder}/{subfamily_filename}.txt'
    internal_id_df = pd.read_csv(internal_id_tbl, sep='\t', index_col= 0)
    old_te_age_folder = f'/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/age_redo'
    old_te_age_tbl = f'{old_te_age_folder}/{subfamily_filename}.txt'
    old_te_age_df = pd.read_csv(old_te_age_tbl, sep='\t')
    te_age_folder = config.te_age_folder
    te_age_tbl = f'{te_age_folder}/{subfamily_filename}.txt'
    te_age_df = pd.read_csv(te_age_tbl, sep='\t')
    te_div_folder = config.kimura_distance_folder
    te_div_tbl = f'{te_div_folder}/{subfamily_filename}.txt'
    te_div_df = pd.read_csv(te_div_tbl, sep='\t')
    subfam_df = filtered_table[filtered_table.repName == subfamily]
    subfam_w_internal_id=pd.merge(subfam_df, internal_id_df, left_index = True, right_on = 'rmsk_index')
    subfam_w_internal_id = subfam_w_internal_id[['genoName','genoStart','genoEnd','strand','repName','repClass','internal_id']]

    subfam_w_internal_id_w_age =pd.merge(subfam_w_internal_id, old_te_age_df, on = 'internal_id')
    subfam_w_internal_id_w_age.rename(columns={'te_age': 'old_te_age'}, inplace=True)
    subfam_w_internal_id_w_age_evi =pd.merge(subfam_w_internal_id_w_age, te_age_df, on = 'internal_id')
    subfam_w_internal_id_w_age_w_div =pd.merge(subfam_w_internal_id_w_age_evi, te_div_df, on = 'internal_id')
    output_folder = '/home/pc575/rds/rds-mi339-kzfps /users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/combined_age_div_evi'
    output_filepath = f'{output_folder}/{subfamily_filename}.txt'
    subfam_w_internal_id_w_age_w_div.to_csv(output_filepath, sep='\t', index = False)
    print(f'done: {subfamily}')
    return subfam_w_internal_id_w_age_w_div
# %%
def main():
    subfamily_list = config.subfamily_list
    with ProcessPoolExecutor(max_workers=40) as executor:
        results = executor.map(make_te_age_bed,subfamily_list)
    combined_te_age=pd.concat(results)
    combined_te_age.to_csv(f'/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/combined_age_div_evi/combined_te_age_div.txt', sep='\t', index=False)
#%%
if __name__ == '__main__':
    main()

# %%
