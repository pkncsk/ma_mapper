#%%
import pyBigWig
import numpy as np
import pandas as pd
import os
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat 
from . import logger
#%%
def extract_bigwig(bigwig_file, chrom, start_list, end_list, strand):
    if isinstance(bigwig_file, str):
        bigwig = pyBigWig.open(bigwig_file)
    elif str(type(bigwig_file)) == "<class 'pyBigWig.bigWigFile'>":
        bigwig = bigwig_file

    bigwig_arrays = []
    for i in range(len(start_list)):
        # zero-based to one-based conversion [BED -> VCF]
        start = start_list[i]
        end = end_list[i]
        bigwig_array = bigwig.values(chrom, start, end, numpy = True)
        if strand == '-':
            bigwig_array = np.flip(bigwig_array)
        bigwig_arrays.append(bigwig_array)
    if strand == '-':
        bigwig_arrays.reverse()
    bigwig_out = np.concatenate(bigwig_arrays)
    return bigwig_out
#%%
def bigwig_io(metadata, bigwig, output_dir = None, save_to_file = False, custom_id = False):
    if isinstance(metadata, str):
        if (os.path.isfile(metadata) == True):
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
    with pyBigWig.open(bigwig) as bigwig_file:
        bigwig_out = []
        for i  in range(len(chrom_list)):
            bigwig_out.append(extract_bigwig(bigwig_file, chrom_list[i], start_list[i], end_list[i], strand_list[i]))

    if save_to_file == True:
        import compress_pickle
        output_filepath = f'{output_dir}/bigwig_out.lzma'
        compress_pickle.dump(bigwig_out, output_filepath, compression="lzma")
        logger.info('done, saving bigwig_out at: ',output_filepath)
    else:
        logger.info('done, returning bigwig_out as object')
        return bigwig_out