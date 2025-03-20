#%%
import sys
sys.path.append('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/ma_mapper/')
from ma_mapper import extract_maf
from ma_mapper import mapper
import pandas as pd
import os
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat 
#%%

#%%
def get_maf_filepath(maf_dir, chrom):
    files = os.listdir(maf_dir)
    maf_filename = f"{chrom}.maf" 
    maf_files = [f for f in files if f.endswith(maf_filename)]

    # Determine the appropriate file to use
    maf_filepath = None
    if any(f.endswith('.maf') for f in maf_files):
        maf_filepath = f"{maf_dir}/{maf_filename}" #.maf is more optimal performance wise 
    elif any(f.endswith('.maf.gz') for f in maf_files):
        maf_filepath = f"{maf_dir}/{maf_filename}.gz" 
    else:
        raise FileNotFoundError(f"No .maf or .maf.gz file found for chromosome {chrom} in {maf_dir}")

    return maf_filepath
#%%
def maf_io(coordinate_table: pd.DataFrame|str, 
           maf:str,
           separated_maf:bool = False, 
           target_species:str = 'Homo_sapiens', 
           count_arg = 'common', 
           save_to_file:bool = False, 
           generate_new_id:bool = False, 
           species_list:list=None,
           e_value_table: str|pd.DataFrame|None=None, 
           internal_id_table: str|pd.DataFrame|None=None, **kwargs):
    if isinstance(coordinate_table, str):
        if os.path.isfile(coordinate_table):
            coordinate_local = pd.read_csv(coordinate_table, sep='\t', header=None)
        else:
            print('coordinate_table file not found')
    else:
        coordinate_local = coordinate_table
    

    print(f'extract from maf target: {maf}')
    if generate_new_id == True:
        meta_id = [f'entry_{index}' for index in coordinate_local.index.astype(str)]
        coordinate_local['meta_id'] = meta_id
    else:
        coordinate_local['meta_id'] = coordinate_local.iloc[:,3]
        meta_id = coordinate_local.meta_id.unique()
    print(meta_id[0])
    grouped = coordinate_local.groupby('meta_id', sort=False)
    chrom_list = grouped.apply(lambda x: x.iloc[:,0].unique()[0]).tolist()
    start_list = grouped.apply(lambda x: x.iloc[:,1].tolist()).tolist()
    end_list = grouped.apply(lambda x: x.iloc[:,2].tolist()).tolist()
    strand_list = grouped.apply(lambda x: x.iloc[:,5].unique()[0]).tolist()
    maf_call_list = []
    for chrom in chrom_list:
        if separated_maf == True:
            print('separated_maf', maf)
            maf_file = get_maf_filepath(maf, chrom)
        else:
            maf_file = maf
        maf_call_list.append(maf_file)
        print(maf_file)
    e_value_df = None
    if isinstance(e_value_table, str):
        e_value_df = pd.read_csv(e_value_table, sep='\t')
    else:
        e_value_df = e_value_table
    internal_id_df = None
    if isinstance(internal_id_table, str):
        internal_id_df = pd.read_csv(internal_id_table, sep='\t')
    else:
        internal_id_df = internal_id_table
    
    with ProcessPoolExecutor(max_workers=40) as executor:
        results = executor.map(extract_maf.extract_maf, meta_id, maf_call_list, chrom_list, start_list, end_list, strand_list, repeat(target_species), repeat(count_arg), repeat(e_value_df), repeat(internal_id_df), repeat(species_list))

    maf_out = []
    for result in results:
        maf_out.append(result)
    if save_to_file == True:
        if isinstance(save_to_file, str):
            output_filepath = save_to_file
        else:
            if isinstance(coordinate_table, str):
                output_dir = '/'.join(str.split(coordinate_table, sep ='/')[:-1])
            else: 
                output_dir = os.path.dirname(os.path.abspath(__file__))
            output_filepath = f'{output_dir}/maf_output.p'
        import compress_pickle
        compress_pickle.dump(maf_out, output_filepath, compression="lzma")
        print('done, saving maf_out at: ', output_filepath)
    else:
        print('done, returning maf_out as object')
        return maf_out

extract_maf.maf_io = maf_io
extract_maf.get_maf_filepath = get_maf_filepath
#%%
alignment_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/alignment/THE1C.fasta.aligned'
alignment_matrix, alignment_coordinate, filters  = mapper.parse_and_filter(alignment_file=alignment_filepath, preprocess_out=True)
#%%
MAF_dir = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/'
maf_matrix = extract_maf.maf_io(coordinate_table=alignment_coordinate, maf = MAF_dir, separated_maf=True, count_arg='common', target_species='hg38')
#%%
coordinate_table = alignment_coordinate
if isinstance(coordinate_table, str):
    if os.path.isfile(coordinate_table):
        coordinate_local = pd.read_csv(coordinate_table, sep='\t', header=None)
    else:
        print('coordinate_table file not found')
else:
    coordinate_local = coordinate_table

maf=MAF_dir
separated_maf=True
generate_new_id=False
e_value_table=None
internal_id_table=None
species_list=None
target_species='hg38'
count_arg='common'
print(f'extract from maf target: {maf}')
if generate_new_id == True:
    meta_id = [f'entry_{index}' for index in coordinate_local.index.astype(str)]
    coordinate_local['meta_id'] = meta_id
else:
    coordinate_local['meta_id'] = coordinate_local.iloc[:,3]
    meta_id = coordinate_local.meta_id.unique()
grouped = coordinate_local.groupby('meta_id', sort=False)
chrom_list = grouped.apply(lambda x: x.iloc[:,0].unique()[0]).tolist()
start_list = grouped.apply(lambda x: x.iloc[:,1].tolist()).tolist()
end_list = grouped.apply(lambda x: x.iloc[:,2].tolist()).tolist()
strand_list = grouped.apply(lambda x: x.iloc[:,5].unique()[0]).tolist()
maf_call_list = []
for chrom in chrom_list:
    if separated_maf == True:
        maf_file = get_maf_filepath(maf, chrom)
    else:
        maf_file = maf
    maf_call_list.append(maf_file)
e_value_df = None
if isinstance(e_value_table, str):
    e_value_df = pd.read_csv(e_value_table, sep='\t')
else:
    e_value_df = e_value_table
internal_id_df = None
if isinstance(internal_id_table, str):
    internal_id_df = pd.read_csv(internal_id_table, sep='\t')
else:
    internal_id_df = internal_id_table

#%%
extract_maf.extract_maf(meta_id[0], maf_call_list[0], chrom_list[0], start_list[0], end_list[0], strand_list[0], target_species, count_arg, e_value_df, internal_id_df, species_list)
#%%
e_value_table_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/e_value/THE1C.txt'
internal_id_table_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/internal_id/THE1C.txt'
maf_matrix_filtered = extract_maf.maf_io(coordinate_table=alignment_coordinate, maf = MAF_dir, separated_maf=True, count_arg='common', e_value_table=e_value_table_filepath, internal_id_table=internal_id_table_filepath)
# %%
