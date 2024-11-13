#%%
import pandas as pd
import numpy as np
import os
import sys
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
from concurrent.futures import ProcessPoolExecutor
#%%
#%%def make_te_age_bed(subfamily):
#subfamily = 'THE1C'
def te_tag(subfamily):
    subfamily_filename = subfamily.replace('/','_')
    print(f'process: {subfamily}')
    #subfamily = 'MER5C'
    #coord_internal_id_folder = config.coord_internal_id_folder
    #te_coord_tbl = f'{coord_internal_id_folder}/{subfamily_filename}.txt'
    #te_coord_df = pd.read_csv(te_coord_tbl, sep='\t')
    internal_id_folder = config.internal_id_folder
    internal_id_tbl = f'{internal_id_folder}/{subfamily_filename}.txt'
    internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
    te_age_folder = config.te_age_folder
    te_age_tbl = f'{te_age_folder}/{subfamily_filename}.txt'
    te_age_df = pd.read_csv(te_age_tbl, sep='\t')
    unique_internal_id = internal_id_df.internal_id.unique()
    age_internal_id=te_age_df.internal_id.unique()
    humaninsertion_folder=config.te_age_human_insertion_folder
    humaninsertion_tbl = f'{humaninsertion_folder}/{subfamily_filename}.txt'
    if os.path.getsize(humaninsertion_tbl) > 1:
        humaninsertion_df = pd.read_csv(humaninsertion_tbl, sep='\t')
        humaninsertion_id=humaninsertion_df.internal_id.unique()
    else:
        print(f"The file {humaninsertion_tbl} is empty.")
        df = pd.DataFrame()  
        humaninsertion_id = []

    segmentaldupe_folder=config.te_age_segmental_folder 
    segmentaldupe_tbl = f'{segmentaldupe_folder}/{subfamily_filename}.txt'
    if os.path.getsize(segmentaldupe_tbl) > 1:
        segmentaldupe_df = pd.read_csv(segmentaldupe_tbl, sep='\t')
        segmentaldupe_id=segmentaldupe_df.internal_id.unique()
    else:
        print(f"The file {segmentaldupe_tbl} is empty.")
        df = pd.DataFrame()  
        segmentaldupe_id = []

    unclassified_folder=config.te_age_unclassified 
    unclassified_tbl = f'{unclassified_folder}/{subfamily_filename}.txt'
    if os.path.getsize(unclassified_tbl) > 1:
        unclassified_df = pd.read_csv(unclassified_tbl, sep='\t')
        unclassified_id=unclassified_df.internal_id.unique()
    else:
        print(f"The file {unclassified_tbl} is empty.")
        df = pd.DataFrame()  
        unclassified_id = []
    tags = []
    for internal_id in unique_internal_id:
        if internal_id not in age_internal_id:
            #print(f'{row.internal_id} is missing')
            tags.append('low_confidence')
        elif internal_id in humaninsertion_id:
            tags.append('human_specific_insertion')
        elif internal_id in segmentaldupe_id:
            tags.append('segmental_duplication')
        elif internal_id in unclassified_id:
            tags.append('unclassified')
        else:
            tags.append('pass')
    dict_prep = {'internal_id': unique_internal_id, 'tag':tags,}
    output_table=pd.DataFrame(dict_prep)
    te_tag_folder = config.te_tag_folder
    output_filepath = f'{te_tag_folder}/{subfamily_filename}.txt'
    output_table.to_csv(output_filepath, sep='\t', index=False)
    print(f'done: {subfamily}')
# %%
def main():
    #for subfamily in ['THE1C']:
    for subfamily in config.subfamily_list:
        te_tag(subfamily)
# %%
if __name__ == '__main__':
    main()
#%%
