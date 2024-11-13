#%%
from unittest import result
import pandas as pd
import numpy as np
import sys
import os
import re
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import itertools
from concurrent.futures import ProcessPoolExecutor
import config_hg38 as config
#%%
def string_to_list(s):
    try:
        return [int(x) for x in re.findall(r'\d+', s)]
    except (ValueError, TypeError):
        print(f"Warning: Could not parse {s}")
        return None

def filter_e(e_table_by_id,e_cutoff = 1e-3):
    age = 0
    segmental = pd.DataFrame()
    unclassified = pd.DataFrame()
    nosig_match = pd.DataFrame()
    internal_id = e_table_by_id.internal_id.unique()[0]
    if e_table_by_id.shape[0] == 1:
        age = np.nan
        e_table_by_id.columns = ['chr_code','divergence','%iden','%gap','BLAST_score','E_value','%iden_flanks','match_total_flank','front_back_matches','both_no_match', 'internal_id']
        e_table_by_id.loc[:,'match_total_flank'] = e_table_by_id.match_total_flank.astype(str).apply(string_to_list)
        if e_table_by_id['match_total_flank'].apply(lambda x: x[0] <= 1).values[0]:
            if e_table_by_id['match_total_flank'].apply(lambda x: x[0] == 1).values[0]:
                segmental = e_table_by_id
            else:
                unclassified = e_table_by_id
            #print(idx, 'possible segmental duplication')
        else:
            nosig_match = e_table_by_id
            #print(idx, 'maybe human specific')
    else:
        
        cutoff_pass_tbl=e_table_by_id[e_table_by_id.E_value.astype('float64')<=e_cutoff].copy()
        #Convert the strings in 'E_val_flanks' to lists
        cutoff_pass_tbl['E_val_flanks'] = cutoff_pass_tbl.E_val_flanks.astype(str).apply(string_to_list)
        # Convert the elements of the lists to floats and replace 'inf' with a large number
        cutoff_pass_tbl['E_val_flanks'] = cutoff_pass_tbl.E_val_flanks.apply(lambda x: [float(i) if i != 'inf' else float('inf') for i in x])
        # Define your threshold
        # Filter rows where at least one element in the list is less than the threshold
        second_pass_tbl = cutoff_pass_tbl[cutoff_pass_tbl.E_val_flanks.apply(lambda x: any(i <= e_cutoff for i in x))]
        te_age=second_pass_tbl.divergence.max()
        age = te_age
    return internal_id, age, segmental, unclassified, nosig_match

def filter_e_for_age(subfamily, e_cutoff = 1e-3):
    e_table = pd.read_csv(f'{config.e_value_folder}/{subfamily}.txt',sep = '\t', low_memory=False)
    output_filepath = f'{config.te_age_folder}/{subfamily}.txt'
    if (os.path.isfile(output_filepath) == False):
        print('start',subfamily)
        grouped =e_table.groupby('internal_id', sort=False)

# Create a list to store the smaller DataFrames
        e_val_table_by_id = [group for _, group in grouped]
       
        
        with ProcessPoolExecutor(max_workers=40) as executor:
            results = executor.map(filter_e, e_val_table_by_id ,itertools.repeat(e_cutoff))                #print(idx, te_age)

        id_list = []
        age_list = []
        nosig_match = []
        segmental = []
        unclassified = []
        for result in results:
            id_list.append(result[0])
            age_list.append(result[1])
            segmental.append(result[2])
            unclassified.append(result[3])
            nosig_match.append(result[4])
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
            unclassified_filepath=f'{config.te_age_unclassified}/{subfamily}.txt'
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

