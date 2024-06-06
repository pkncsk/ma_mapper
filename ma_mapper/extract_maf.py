#%%
import numpy as np
import pandas as pd
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat 
from Bio.AlignIO import MafIO
from . import logger
if sys.version_info >= (3, 8, 0):
    from typing import Literal, Tuple, List
else:
    from typing_extensions import Literal, Tuple, List
#%%
_AGEARG = Literal[None,'extract','calibrate']
_COUNTARG = Literal['human_ref','coverage','common', 'common_raw', 'total_raw','a','t','c','g','base_freq']
def extract_maf(maf_file:str, 
                chrom:str, 
                start:int, 
                end:int, 
                strand:str,
                target_species:str = 'Homo_sapiens', 
                count_arg:_COUNTARG = 'common',
                age:int=0, 
                age_arg:_AGEARG = None, 
                species_age_table_file:pd.DataFrame = None):
    #print(target_species,chrom, start, end, strand)
    if age_arg is not None:
        if isinstance(species_age_table_file, str):
            if (os.path.isfile(species_age_table_file) == True):
                species = pd.read_table(species_age_table_file, sep='\t')
            else:
                logger.error('metadata file not found')
        else:
            species = species_age_table_file
    maf_id = f'{target_species}.{chrom}'

    index_maf = MafIO.MafIndex(f'{maf_file}.mafindex', maf_file, maf_id) 
    n_strand = -1 if strand == '-' else 1
    results =index_maf.get_spliced(start,end,n_strand)
    if age_arg is None:
        collector = {seqrec.id.split('.')[0]: np.array(list(seqrec.seq.lower())) for seqrec in results}
    else:
        if age_arg == 'calibrate':
            collector = {seqrec.id.split('.')[0]: np.array(list(seqrec.seq.lower())) for seqrec in results if age >= species[species.meta_name == seqrec.id.split('.')[0]]['Estimated Time (MYA)'].to_list()[0]}
        elif age_arg == 'extract':
            for seqrec in results:
                species_age = species[species.meta_name == seqrec.id.split('.')[0]]['Estimated Time (MYA)'].to_list()[0]
                if species_age >= age:
                    age = species_age
            return age

    
    array_transposed=np.array(list(collector.values())).transpose()
    output_array=[]
    ref_alleles=collector['Homo_sapiens']
    for idx,ref_pos in enumerate(array_transposed):
        (unique, counts) = np.unique(ref_pos, return_counts=True)
        frequencies = dict(zip(unique, counts))
        ref_allele=ref_alleles[idx]
        #frequencies.pop('-', None) #->count deletion/insertion
        total = sum(frequencies.values())
        if count_arg == 'human_ref':
            alt_count = total - frequencies[ref_allele]
            alt_freq=alt_count/total
            output_array.append(alt_freq)
        elif count_arg in ['coverage', 'total_raw']:
            output_array.append(total)
        elif count_arg in ['common', 'common_raw']:
            common_allele = max(frequencies, key=frequencies.get)
            common_count = frequencies.get(common_allele)
            if count_arg == 'common':
                common_freq = common_count/total
                output_array.append(common_freq)
            else:
                output_array.append(common_count)
        elif count_arg in ['a','t','c','g']:
            nucl_count = frequencies.get(count_arg, 0)
            output_array.append(nucl_count)
        elif count_arg == 'base_freq':
            output_array.append(frequencies)
    return output_array
#%%
def maf_io(metadata: pd.DataFrame|str, 
           maf:str,
           output_dir:str|None = None, 
           separated_maf:bool = False, 
           target_species:str = 'Homo_sapiens', 
           count_arg:_COUNTARG = 'common', 
           save_to_file:bool = False, 
           custom_id:bool = False, 
           age_arg:_AGEARG = None, 
           species_age_table_file:str|None = None, 
           age_list: List|str|pd.DataFrame|None=None, **kwargs):
    if isinstance(metadata, str):
        if os.path.isfile(metadata):
            metadata_local = pd.read_csv(metadata, sep='\t')
        else:
            logger.error('metadata file not found')
    else:
        metadata_local = metadata
    if output_dir is None:
        if isinstance(metadata, str):
            output_dir = '/'.join(str.split(metadata, sep ='/')[:-1])
        else: 
            output_dir = os.path.dirname(os.path.abspath(__file__))

    logger.info(f'extract from maf target: {maf}')
    if custom_id == False:
        meta_id = [f'sample_n{index}' for index in metadata_local.index.astype(str)]
        metadata_local['meta_id'] = meta_id
    else:
        metadata_local['meta_id'] = metadata_local.iloc[:,4]
        meta_id = metadata_local.meta_id.unique()

    grouped = metadata_local.groupby('meta_id', sort=False)
    chrom_list = grouped.apply(lambda x: x.iloc[:,0].unique()[0], include_groups=False).tolist()
    start_list = grouped.apply(lambda x: x.iloc[:,1].tolist(), include_groups=False).tolist()
    end_list = grouped.apply(lambda x: x.iloc[:,2].tolist(), include_groups=False).tolist()
    strand_list = grouped.apply(lambda x: x.iloc[:,3].unique()[0], include_groups=False).tolist()
    maf_call_list = []
    for chrom in chrom_list:
        if separated_maf == True:
            maf_file = f'{maf}.{chrom}'
        else:
            maf_file = maf
        maf_call_list.append(maf_file)

    if age_arg is None:
        logger.info('normal maf extraction without calibration')
        with ProcessPoolExecutor(max_workers=40) as executor:
            results  = executor.map(extract_maf, maf_call_list, chrom_list, start_list, end_list, strand_list, repeat(target_species) ,repeat(count_arg))
    elif age_arg == 'extract':
        logger.info('extract age depth for further use')
        with ProcessPoolExecutor(max_workers=40) as executor:
            results  = executor.map(extract_maf, maf_call_list, chrom_list, start_list, end_list, strand_list, repeat(target_species) ,repeat(count_arg), repeat(age_list), repeat(age_arg), repeat(species_age_table_file))
    elif age_arg == 'calibrate':
        if age_list is None:
            logger.info('calibrate age: no age_list provided, fetching from metadata')
            age_list_local = []
            if 'te_age' in metadata_local:
                for uniq_meta_id in meta_id:
                    metadata_by_id = metadata_local[metadata_local.meta_id ==uniq_meta_id]
                    age_list_local.append(metadata_by_id.te_age.unique()[0])
            else:
                logger.info('no te_age infomation from metadata, running match_age_to_id_metadata()')
                from . import mapper
                metadata_age=mapper.match_age_to_id_metadata(metadata_local, **kwargs)
                age_list_local=metadata_age.drop_duplicates(subset='meta_id').te_age.to_list()
        else:
            logger.info('calibrate age: using preexisting age_list')
            if isinstance(age_list, List):
                age_list_local = age_list
            elif isinstance(age_list, str):
                age_list_local_df = pd.read_csv(age_list, sep='\t')
                age_list_local = age_list_local_df['te_age']
            elif isinstance(age_list, pd.DataFrame):
                age_list_local = age_list['te_age'].to_list()
        with ProcessPoolExecutor(max_workers=40) as executor:
            results  = executor.map(extract_maf, maf_call_list, chrom_list, start_list, end_list, strand_list, repeat(target_species) ,repeat(count_arg), age_list_local, repeat(age_arg), repeat(species_age_table_file))
    maf_out = []
    for result in results:
        maf_out.append(result)
    if save_to_file == True:
        import compress_pickle
        output_filepath = f'{output_dir}/maf_out.lzma'
        compress_pickle.dump(maf_out, output_filepath, compression="lzma")
        logger.info('done, saving maf_out at: ', output_filepath)
    else:
        logger.info('done, returning maf_out as object')
        return maf_out
#%%
