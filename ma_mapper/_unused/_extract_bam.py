#%%
import math
import numpy as np
import pandas as pd
import pysam
import os
from .. import logger
#-> BAM file overlay 
import sys
if sys.version_info >= (3, 8, 0):
    from typing import Literal, Tuple, List
else:
    from typing_extensions import Literal, Tuple, List
def normal_array(width=1, sigma=1, odd=0):
    ''' Returns an array of the normal distribution of the specified width '''
    sigma2 = sigma ** 2
    def normal_func(x):
        return math.exp(-x * x / (2 * sigma2))
    values = [normal_func(x) for x in range(-width, width + odd)]
    return np.array(values)

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

        if bam_format == 'bam_min':
            output = np.minimum(profile_normals_forward, profile_normals_reverse)
        elif bam_format == 'bam_max':
            output = np.maximum(profile_normals_forward, profile_normals_reverse)
        elif bam_format == 'bam_forward':
            output = profile_reads_forward
        elif bam_format == 'bam_reverse':
            output = profile_reads_reverse

        if strand == '-':
            output = np.flip(output)

        output = output[probe_length:probe_length+fragment_length]
        result_array.append(output)

    np_result = np.concatenate(result_array) if result_array else np.zeros(np.sum(end_list - start_list), dtype=np.uint8)
    return np_result

def bam_io(metadata: str|pd.DataFrame, 
              bam_file:str,
              bam_format:_FORMAT, 
              output_dir:str = None, 
              save_to_file:bool = False, 
              custom_id:bool = False, 
              **kwargs)->List:
    if isinstance(metadata, str):
        if (os.path.isfile(metadata) == True):
            metadata = pd.read_csv(metadata, sep='\t')
        else:
            print('metadata file not found')
    else:
        metadata = metadata
    if output_dir is None:
        if isinstance(metadata, str):
            output_dir = '/'.join(str.split(metadata, sep ='/')[:-1])
        else:
            output_dir = os.path.dirname(os.path.abspath(__file__))
    if custom_id == False:
        meta_id = 'sample_n'+ metadata.index.astype(str)
        metadata['meta_id'] = meta_id
    else:
        metadata['meta_id'] = metadata.id
        meta_id = metadata.id.unique()   
    grouped = metadata.groupby('meta_id')

    chrom_list = grouped.apply(lambda x: x.iloc[:,0].unique()[0], include_groups=False).tolist()
    start_list = grouped.apply(lambda x: x.iloc[:,1].tolist(), include_groups=False).tolist()
    end_list = grouped.apply(lambda x: x.iloc[:,2].tolist(), include_groups=False).tolist()
    strand_list = grouped.apply(lambda x: x.iloc[:,3].unique()[0], include_groups=False).tolist()
    eb_kwargs = {key: value for key, value in kwargs.items() if key in ('offset','probe_length', 'smoothing_length')}
    bam = []
    with pysam.AlignmentFile(bam_file, "rb") as bam_file:
        for i  in range(len(chrom_list)):
            np_result = extract_bam(bam_file,chrom_list[i], start_list[i], end_list[i], strand_list[i], bam_format, **eb_kwargs)
            bam.append(np_result)

    import compress_pickle
    if save_to_file == True:
        output_filepath = output_dir+'/bam_'+bam_format+'.lzma'
        compress_pickle.dump(bam, output_filepath, compression="lzma")
        logger.info('done, saving bam_min at: '+output_filepath)

    else:
        logger.info('done, returning mapped bam as a list. format: %s', bam_format)
        return bam