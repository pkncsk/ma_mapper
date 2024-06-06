#%%
import pandas as pd
import sys
import os
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/stable')
import config_hg38_repeatlib as config
from concurrent.futures import ProcessPoolExecutor
#%%
#subfamily ='THE1C'
def make_te_gap(subfamily):
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
    subfam_df = filtered_table[filtered_table.repName == subfamily]
    subfam_w_internal_id=pd.merge(subfam_df, internal_id_df, left_index = True, right_on = 'rmsk_index')
    #subfam_w_internal_id = subfam_w_internal_id[['genoName','genoStart','genoEnd','strand','repName','repClass','internal_id']]
    subfam_w_internal_id['gap'] = subfam_w_internal_id.groupby('internal_id')['genoStart'].shift(-1) - subfam_w_internal_id['genoEnd']
    gap_tbl=subfam_w_internal_id[['internal_id','gap']]
    
    #gap_tbl[~gap_tbl.gap.isna()]
    output_folder = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/gap'
    output_filepath = f'{output_folder}/{subfamily_filename}.txt'
    gap_tbl.to_csv(output_filepath, sep='\t', index = False)
    gap_tbl['internal_id'] = subfamily+'_'+gap_tbl.internal_id.astype(str)
    return gap_tbl
# %%
def main():
    subfamily_list = config.subfamily_list
    with ProcessPoolExecutor(max_workers=40) as executor:
        results = executor.map(make_te_gap,subfamily_list)
    combined_te_gap=pd.concat(results)
    combined_te_gap.to_csv(f'/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/gap/combined_te_gap.txt', sep='\t', index=False)
#%%
if __name__ == '__main__':
    main()
#%%