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
        #print(start,end)
        bigwig_array = bigwig.values(chrom, start, end, numpy = True)
        if strand == '-':
            bigwig_array = np.flip(bigwig_array)
        bigwig_arrays.append(bigwig_array)
    if strand == '-':
        bigwig_arrays.reverse()
    bigwig_out = np.concatenate(bigwig_arrays)
    return bigwig_out
#%%
def bigwig_io(coordinate_table, 
              bigwig, 
              output_filepath = None, 
              save_to_file = False, 
              custom_id = False,
              custom_prefix='entry'):
    if isinstance(coordinate_table, str):
        if (os.path.isfile(coordinate_table) == True):
            coordinate_local = pd.read_csv(coordinate_table, sep='\t', header=None)
        else:
            logger.error('coordinate_table file not found')
    else:
        coordinate_local = coordinate_table
    if output_filepath is None:
        if isinstance(coordinate_table, str):
            output_filepath = '/'.join(str.split(coordinate_table, sep ='/')[:-1])
        else:
            output_filepath = os.path.dirname(os.path.abspath(__file__))
    if custom_id == False:
        meta_id = [f'{custom_prefix}_{index}' for index in coordinate_local.index.astype(str)]
        coordinate_local['meta_id'] = meta_id
    else:
        coordinate_local['meta_id'] = coordinate_local.iloc[:,3]
        meta_id = coordinate_local.meta_id.unique()    
    grouped = coordinate_local.groupby('meta_id', sort=False)
    chrom_list = grouped.apply(lambda x: x.iloc[:,0].unique()[0], include_groups=False).tolist()
    start_list = grouped.apply(lambda x: x.iloc[:,1].tolist(), include_groups=False).tolist()
    end_list = grouped.apply(lambda x: x.iloc[:,2].tolist(), include_groups=False).tolist()
    strand_list = grouped.apply(lambda x: x.iloc[:,5].unique()[0], include_groups=False).tolist()
    with pyBigWig.open(bigwig) as bigwig_file:
        bigwig_out = []
        for i  in range(len(chrom_list)):
            bigwig_out.append(extract_bigwig(bigwig_file, chrom_list[i], start_list[i], end_list[i], strand_list[i]))

    if save_to_file == True:
        import compress_pickle
        compress_pickle.dump(bigwig_out, output_filepath, compression="lzma")
        logger.info('done, saving bigwig_out at: ',output_filepath)
    else:
        logger.info('done, returning bigwig_out as object')
        return bigwig_out