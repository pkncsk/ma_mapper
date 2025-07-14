#%%
import math
import numpy as np
import pandas as pd
import pysam
import os
from . import logger
#-> BAM file overlay 
import sys
if sys.version_info >= (3, 8, 0):
    from typing import Literal, Tuple, List
else:
    from typing_extensions import Literal, Tuple, List
def normal_array(width=1, sigma=1, odd=0):
#''' Returns an array of the normal distribution of the specified width '''
    sigma2 = sigma ** 2
    def normal_func(x):
        return math.exp(-x * x / (2 * sigma2))
    # width is the half of the distribution
    values = list(map(normal_func, range(-width, width + odd)))
    normal_array = np.array(values)
    return normal_array
#%%
_FORMAT = Literal['read_max','read_min','read_forward','read_reverse', 'normal_max','normal_min','normal_forward','normal_reverse','read_sum']

def extract_bam(bam_file:str, chrom:str, start_list:List, end_list:List, strand:str, bam_format:_FORMAT, offset:int = 5, probe_length:int =100, smoothing_length:int= 100, drop_dup=False, drop_unmap=False, drop_2=False, drop_multimap=False)->np.ndarray:
    if isinstance(bam_file, str):
        bam_file = pysam.AlignmentFile(bam_file, "rb")
    normal = normal_array(width=probe_length, sigma=smoothing_length, odd=1)

    result_array = []
    for i in range(len(start_list)):
        start = start_list[i]
        end = end_list[i]
        fragment_length = end - start
        profile_normals_forward = np.zeros(fragment_length + (2*probe_length), dtype=np.float16)
        profile_normals_reverse = np.zeros(fragment_length + (2*probe_length), dtype=np.float16)
        profile_reads_forward = np.zeros(fragment_length + (2*probe_length), dtype=np.uint8)
        profile_reads_reverse = np.zeros(fragment_length + (2*probe_length), dtype=np.uint8)
        mapped_reads = bam_file.fetch(chrom, start, end)
        
        for mapped_read in mapped_reads:
            if drop_dup and mapped_read.is_duplicate:
                continue
            if drop_unmap and mapped_read.is_unmapped:
                continue
            if drop_2 and (mapped_read.is_secondary or mapped_read.is_supplementary):
                continue
            if drop_multimap and mapped_read.has_tag("NH") and mapped_read.get_tag("NH") > 1:
                continue
            if mapped_read.is_reverse is False:
                read_position = mapped_read.reference_start - start + offset + probe_length
            else:
                read_position = mapped_read.reference_end - start - offset - 1 + probe_length

            start_in_window = read_position - probe_length
            end_in_window = read_position + probe_length + 1

            if (start_in_window < 0) or (end_in_window > fragment_length + (2*probe_length)):
                continue
            
            if strand == '+':
                if mapped_read.is_reverse is False:
                    np.add(profile_normals_forward[start_in_window:end_in_window], normal, out=profile_normals_forward[start_in_window:end_in_window])
                    profile_reads_forward[read_position] += 1
                else:
                    np.add(profile_normals_reverse[start_in_window:end_in_window], normal, out=profile_normals_reverse[start_in_window:end_in_window])
                    profile_reads_reverse[read_position] += 1
            else:
                if mapped_read.is_reverse is False:
                    np.add(profile_normals_reverse[start_in_window:end_in_window], normal, out=profile_normals_reverse[start_in_window:end_in_window])
                    profile_reads_reverse[read_position] += 1
                else:
                    np.add(profile_normals_forward[start_in_window:end_in_window], normal, out=profile_normals_forward[start_in_window:end_in_window])
                    profile_reads_forward[read_position] += 1

        if bam_format == 'normal_min':
            output = np.minimum(profile_normals_forward, profile_normals_reverse)
        elif bam_format == 'normal_max':
            output = np.maximum(profile_normals_forward, profile_normals_reverse)
        elif bam_format == 'read_min':
            output = np.minimum(profile_reads_forward, profile_reads_reverse)
        elif bam_format == 'read_max':
            output = np.maximum(profile_reads_forward, profile_reads_reverse)
        elif bam_format == 'normal_forward':
            output = np.array(profile_normals_forward)
        elif bam_format == 'normal_reverse':
            output = np.array(profile_normals_reverse)
        elif bam_format == 'read_forward':
            output = np.array(profile_reads_forward)
        elif bam_format == 'read_reverse':
            output = np.array(profile_reads_reverse)
        elif bam_format == 'read_sum':
            output = np.add(profile_reads_forward, profile_reads_reverse)
        if strand == '-':
            output = np.flip(output)

        output = output[probe_length:probe_length+fragment_length]
        result_array.append(output)

    np_result = np.concatenate(result_array) if result_array else np.zeros(np.sum(np.subtract(end_list,start_list)), dtype=np.uint8)
    return np_result
#%%
def bam_io(coordinate_table: str|pd.DataFrame, 
              bam_file:str,
              bam_format:_FORMAT, 
              save_to_file:bool = False, 
              generate_new_id:bool = False,
              **kwargs)->List:
    

    if isinstance(coordinate_table, str):
        if (os.path.isfile(coordinate_table) == True):
            metadata_local = pd.read_csv(coordinate_table, sep='\t', header=None)
        else:
            logger.error('coordinate_table file not found')
    else:
        metadata_local = coordinate_table
    if generate_new_id == True:
        meta_id = [f'entry_{index}' for index in metadata_local.index.astype(str)]
        metadata_local['meta_id'] = meta_id
    else:
        metadata_local['meta_id'] = metadata_local.iloc[:,3]
        meta_id = metadata_local.meta_id.unique()

    with pysam.AlignmentFile(bam_file, "rb") as bam_file:
        eb_kwargs = {key: value for key, value in kwargs.items() if key in ('offset','probe_length', 'smoothing_length')}
        grouped = metadata_local.groupby('meta_id', sort=False)
        chrom_list = grouped.apply(lambda x: x.iloc[:,0].unique()[0]).tolist()
        start_list = grouped.apply(lambda x: x.iloc[:,1].tolist()).tolist()
        end_list = grouped.apply(lambda x: x.iloc[:,2].tolist()).tolist()
        strand_list = grouped.apply(lambda x: x.iloc[:,5].unique()[0]).tolist()
        bam = []
        for i  in range(len(chrom_list)):
            np_result = extract_bam(bam_file,chrom_list[i], start_list[i], end_list[i], strand_list[i], bam_format, **eb_kwargs)
            bam.append(np_result)

    if save_to_file:
        if isinstance(save_to_file, str)==True:
            output_filepath = save_to_file
        else:
            if isinstance(coordinate_table, str)==True:
                output_dir = '/'.join(str.split(coordinate_table, sep ='/')[:-1])
            else:
                output_dir = os.path.dirname(os.path.abspath(__file__))
            output_filepath = f'{output_dir}/bam_{bam_format}.p'
        import compress_pickle
        compress_pickle.dump(bam, output_filepath, compression="lzma")
        logger.info('done, saving bam_min at: ',output_filepath)

    else:
        logger.info('done, returning mapped bam as a list. format: %s', bam_format)
        return bam