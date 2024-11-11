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
    values = np.array(values)
    return values
#%%
_FORMAT = Literal['read_max','read_min','read_forward','read_reverse', 'normal_max','normal_min','normal_forward','normal_reverse','read_sum']
def _extract_bam(bam_file:str, 
                chrom:str, 
                start_list:List, 
                end_list:List, 
                strand:str, 
                bam_format:_FORMAT,
                offset:int = 5, 
                probe_length:int =100, 
                smoothing_length:int= 100)->np.ndarray:
    #print(chrom, start_list, end_list, strand)
    if isinstance(bam_file, str):
        bam_file = pysam.AlignmentFile(bam_file, "rb")
    normal = normal_array(width=probe_length, sigma=smoothing_length, odd=1)
    result_array = []
    for i in range(len(start_list)):
        start = start_list[i]
        end = end_list[i]
        fragment_length = end - start
        print(fragment_length)
        profile_normals_forward = np.zeros(fragment_length + (2*probe_length), dtype=np.float16)
        profile_normals_reverse = np.zeros(fragment_length + (2*probe_length), dtype=np.float16)
        profile_reads_forward = np.zeros(fragment_length + (2*probe_length), dtype=np.uint8)
        profile_reads_reverse = np.zeros(fragment_length + (2*probe_length), dtype=np.uint8)
        mapped_reads = bam_file.fetch(chrom, start, end)
        
        for mapped_read in mapped_reads:
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
                    profile_normals_forward[start_in_window:end_in_window] += normal
                    profile_reads_forward[read_position] += 1
                else:
                    profile_normals_reverse[start_in_window:end_in_window] += normal
                    profile_reads_reverse[read_position] += 1
            #for - strand, throw forward read into reverse bin
            elif strand == '-':
                if mapped_read.is_reverse is False:
                    profile_normals_reverse[start_in_window:end_in_window] += normal
                    profile_reads_reverse[read_position] += 1
                else:
                    profile_normals_forward[start_in_window:end_in_window] += normal
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

    #failsafe in case of no mapped reads
    if len(result_array) != 0:
        if strand == '-':
            result_array.reverse()
        np_result = np.concatenate(result_array)
    else:
        np_result = np.zeros(np.sum(np.subtract(end_list,start_list)), dtype=np.uint8)

    return np_result

def extract_bam(bam_file:str, chrom:str, start_list:List, end_list:List, strand:str, bam_format:_FORMAT, offset:int = 5, probe_length:int =100, smoothing_length:int= 100)->np.ndarray:
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
_IO = Literal['legacy','optimized']
_MODE = Literal['legacy','optimized']
def bam_io(metadata: str|pd.DataFrame, 
              bam_file:str,
              bam_format:_FORMAT, 
              output_dir:str = None, 
              save_to_file:bool = False, 
              custom_id:bool = False,
              custom_prefix = 'entry',
              io:_IO = 'optimized',
              mode:_MODE = 'optimized', 
              **kwargs)->List:
    

    if isinstance(metadata, str):
        if (os.path.isfile(metadata) == True):
            metadata_local = pd.read_csv(metadata, sep='\t', header=None)
        else:
            logger.error('metadata file not found')
    else:
        metadata_local = metadata
    if output_dir is None:
        if isinstance(metadata, str):
            output_dir = '/'.join(str.split(metadata, sep ='/')[:-1])
        else:
            output_dir = os.path.dirname(os.path.abspath(__file__))
    if custom_id == False:
        meta_id = [f'{custom_prefix}_{index}' for index in metadata_local.index.astype(str)]
        metadata_local['meta_id'] = meta_id
    else:
        metadata_local['meta_id'] = metadata_local.iloc[:,3]
        meta_id = metadata_local.meta_id.unique()

    with pysam.AlignmentFile(bam_file, "rb") as bam_file:
        eb_kwargs = {key: value for key, value in kwargs.items() if key in ('offset','probe_length', 'smoothing_length')}
        if io == 'optimized': 
            grouped = metadata_local.groupby('meta_id', sort=False)
            chrom_list = grouped.apply(lambda x: x.iloc[:,0].unique()[0], include_groups=False).tolist()
            start_list = grouped.apply(lambda x: x.iloc[:,1].tolist(), include_groups=False).tolist()
            end_list = grouped.apply(lambda x: x.iloc[:,2].tolist(), include_groups=False).tolist()
            strand_list = grouped.apply(lambda x: x.iloc[:,5].unique()[0], include_groups=False).tolist()
            bam = []
            for i  in range(len(chrom_list)):
                if mode == 'legacy':
                    np_result = _extract_bam(bam_file,chrom_list[i], start_list[i], end_list[i], strand_list[i], bam_format, **eb_kwargs)
                elif mode == 'optimized':
                    np_result = extract_bam(bam_file,chrom_list[i], start_list[i], end_list[i], strand_list[i], bam_format, **eb_kwargs)
                bam.append(np_result)
        elif io == 'legacy':
            bam = []
            for uniq_meta_id in meta_id:
                metadata_by_id = metadata_local[metadata_local.meta_id == uniq_meta_id]
                chrom = metadata_by_id.iloc[:,0].unique()[0]
                start_list = metadata_by_id.iloc[:,1].to_list()
                end_list = metadata_by_id.iloc[:,2].to_list()
                strand = metadata_by_id.iloc[:,5].unique()[0]
                if mode == 'legacy':
                    np_result = _extract_bam(bam_file, chrom, start_list, end_list, strand, bam_format, **eb_kwargs)
                elif mode == 'optimized':
                    np_result = extract_bam(bam_file,chrom, start_list, end_list, strand, bam_format, **eb_kwargs)
                bam.append(np_result)

    if save_to_file == True:
        import compress_pickle
        output_filepath = f'{output_dir}/bam_{bam_format}.lzma'
        compress_pickle.dump(bam, output_filepath, compression="lzma")
        logger.info('done, saving bam_min at: ',output_filepath)

    else:
        logger.info('done, returning mapped bam as a list. format: %s', bam_format)
        return bam