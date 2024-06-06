#%%
import pandas as pd
import sys
import os
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/stable')
import config_hg38_repeatlib as config
from concurrent.futures import ProcessPoolExecutor
#%%
def make_te_age_bed(subfamily):
    #subfamily = 'THE1C'
    subfamily_filename = subfamily.replace('/','_')
    output_folder = config.te_age_bed_folder
    output_filepath = f'{output_folder}/{subfamily_filename}.txt'
    print(f'process: {subfamily}')
    #subfamily = 'MER5C'
    #coord_internal_id_folder = config.coord_internal_id_folder
    #te_coord_tbl = f'{coord_internal_id_folder}/{subfamily_filename}.txt'
    #te_coord_df = pd.read_csv(te_coord_tbl, sep='\t')
    filtered_table = config.filtered_table
    internal_id_folder = config.internal_id_folder
    internal_id_tbl = f'{internal_id_folder}/{subfamily_filename}.txt'
    internal_id_df = pd.read_csv(internal_id_tbl, sep='\t', index_col= 0)
    te_age_folder = config.te_age_folder
    te_age_tbl = f'{te_age_folder}/{subfamily_filename}.txt'
    te_age_df = pd.read_csv(te_age_tbl, sep='\t', index_col= 0)
    
    subfam_df = filtered_table[filtered_table.repName == subfamily]
    subfam_w_internal_id=pd.merge(subfam_df, internal_id_df, left_index = True, right_on = 'rmsk_index')
    subfam_w_internal_id = subfam_w_internal_id[['genoName','genoStart','genoEnd','strand','repName','internal_id']]
    
    subfam_w_internal_id_w_age =pd.merge(subfam_w_internal_id, te_age_df, on = 'internal_id')
    
    bed_output=pd.DataFrame(subfam_w_internal_id_w_age[['genoName','genoStart','genoEnd']])
    bed_output.columns = ['chrom','start','end']
    bed_output['name'] = subfam_w_internal_id_w_age.repName + '_' + subfam_w_internal_id_w_age.internal_id.astype(str)
    bed_output['score'] = subfam_w_internal_id_w_age['te_age']
    bed_output['strand'] = subfam_w_internal_id_w_age['strand']
    bed_output = bed_output.fillna(0)
    bed_output.to_csv(f'{config.te_age_bed_folder}/{subfamily}.txt', sep='\t', index=False, header=False)
    print(f'done: {subfamily}')
    return bed_output
# %%
def main():
    subfamily_list = config.subfamily_list
    with ProcessPoolExecutor(max_workers=40) as executor:
        results = executor.map(make_te_age_bed,subfamily_list)
    combined_te_age=pd.concat(results)
    combined_te_age.to_csv(f'{config.te_age_bed_folder}/combined_te_age.txt', sep='\t', index=False, header=False)
#%%
if __name__ == '__main__':
    main()
#%%