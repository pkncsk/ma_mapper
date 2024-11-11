#%%
import pandas as pd
import numpy as np
import sys
import os
import ast
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
#%%
def filter_e_for_age(subfamily, e_cutoff = 1e-3):
    e_table = pd.read_csv(f'{config.e_value_folder}/{subfamily}.txt',sep = '\t', low_memory=False)
    output_filepath = f'{config.te_age_folder}/{subfamily}.txt'
    if (os.path.isfile(output_filepath) == False):
        print('start',subfamily)
        internal_id=e_table.internal_id.unique()
        print(f'total: {len(internal_id)}')
        id_list = []
        age_list = []
        nosig_match = []
        segmental = []
        unclassified = []
        for idx in internal_id:
            print(idx)
            id_list.append(idx)
            e_table_by_id=e_table[e_table.internal_id==idx]
            if e_table_by_id.shape[0] == 1:
                age_list.append(np.nan)
                e_table_by_id.columns = ['chr_code','divergence','%iden','%gap','BLAST_score','E_value','%iden_flanks','match_total_flank','front_back_matches','both_no_match', 'internal_id']
                e_table_by_id.loc[:,'match_total_flank'] = e_table_by_id.match_total_flank.apply(ast.literal_eval)
                if e_table_by_id['match_total_flank'].apply(lambda x: x[0] <= 1).values[0]:
                    if e_table_by_id['match_total_flank'].apply(lambda x: x[0] == 1).values[0]:
                        segmental.append(e_table_by_id)
                    else:
                        unclassified.append(e_table_by_id)
                    #print(idx, 'possible segmental duplication')
                else:
                    nosig_match.append(e_table_by_id)
                    #print(idx, 'maybe human specific')
            else:
                
                cutoff_pass_tbl=e_table_by_id[e_table_by_id.E_value.astype('float64')<=e_cutoff].copy()
                #Convert the strings in 'E_val_flanks' to lists
                cutoff_pass_tbl['E_val_flanks'] = cutoff_pass_tbl.E_val_flanks.apply(ast.literal_eval)
                # Convert the elements of the lists to floats and replace 'inf' with a large number
                cutoff_pass_tbl['E_val_flanks'] = cutoff_pass_tbl.E_val_flanks.apply(lambda x: [float(i) if i != 'inf' else float('inf') for i in x])
                # Define your threshold
                # Filter rows where at least one element in the list is less than the threshold
                second_pass_tbl = cutoff_pass_tbl[cutoff_pass_tbl.E_val_flanks.apply(lambda x: any(i < e_cutoff for i in x))]
                te_age=second_pass_tbl.divergence.max()
                age_list.append(te_age)
                #print(idx, te_age)

        dict_prep = {'internal_id': id_list, 'te_age':age_list,}
        output_table=pd.DataFrame(dict_prep)
        
        output_table.to_csv(output_filepath, sep='\t', index=False)

        if nosig_match:
            nosig_match_df=pd.concat(nosig_match)
            nosig_filepath=f'{config.te_age_human_insertion_folder}/{subfamily}.txt'
            nosig_match_df.to_csv(nosig_filepath, sep='\t', index=False)
        if segmental:
            segmental_df=pd.concat(segmental)
            segmental_filepath=f'{config.te_age_segmental_folder}/{subfamily}.txt'
            segmental_df.to_csv(segmental_filepath, sep='\t', index=False)
        
        if unclassified:
            unclassified_df=pd.concat(unclassified)
            unclassified_filepath=f'{config.te_age_segmental_folder}/{subfamily}.txt'
            unclassified_df.to_csv(unclassified_filepath, sep='\t', index=False)
        
        print('done',subfamily)
    else:
        print('already done', subfamily)
# %%
def main():
    #for subfamily in ['THE1C']:
    for subfamily in config.subfamily_list:
        filter_e_for_age(subfamily)
# %%
if __name__ == '__main__':
    main()
#%%
# %%
