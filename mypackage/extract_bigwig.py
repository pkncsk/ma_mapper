#%%
import pyBigWig
import numpy as np
import pandas as pd
import os
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat 
#%%
def extract_bigwig(bigwig_file, chrom, start_list, end_list, strand):
    bigwig_arrays = []
    for i in range(len(start_list)):
        # zero-based to one-based conversion [BED -> VCF]
        start = start_list[i]
        end = end_list[i]

        bigwig = pyBigWig.open(bigwig_file)
        bigwig_array = bigwig.values(chrom, start, end, numpy = True)
        if strand == '-':
            bigwig_array = np.flip(bigwig_array)
        bigwig_arrays.append(bigwig_array)
    if strand == '-':
        bigwig_arrays.reverse()
    bigwig_out = np.concatenate(bigwig_arrays)
    return bigwig_out
#%%
def fetch_bigwig(metadata_input, bigwig_input, output_dir = None, save_to_file = False, custom_id = False):
    if isinstance(metadata_input, str):
        if (os.path.isfile(metadata_input) == True):
            metadata = pd.read_csv(metadata_input, sep='\t')
        else:
            print('metadata file not found')
    else:
        metadata = metadata_input
    if output_dir is None:
        if isinstance(metadata_input, str):
            output_dir = '/'.join(str.split(metadata_input, sep ='/')[:-1])
        else:
            output_dir = os.path.dirname(os.path.abspath(__file__))
    import logging
    if custom_id == False:
        meta_id = 'sample_n'+ metadata.index.astype(str)
        metadata['meta_id'] = meta_id
    else:
        metadata['meta_id'] = metadata.iloc[:,4]
        meta_id = metadata.iloc[:,4].unique()    
    chrom_list = []
    start_list = []
    end_list = []
    strand_list = []    
    for uniq_meta_id in meta_id:
        metadata_by_id = metadata[metadata.meta_id == uniq_meta_id]
        chrom = metadata_by_id.iloc[:,0].unique()[0]
        chrom_list.append(chrom)
        start_list.append(metadata_by_id.iloc[:,1].to_list())
        end_list.append(metadata_by_id.iloc[:,2].to_list())
        strand_list.append(metadata_by_id.iloc[:,3].unique()[0])
    with ProcessPoolExecutor(max_workers=40) as executor:
        results  = executor.map(extract_bigwig, repeat(bigwig_input), chrom_list, start_list,end_list,strand_list)
    bigwig_out = []
    for result in results:
        bigwig_out.append(result)
    if save_to_file == True:
        import compress_pickle
        output_filepath = output_dir+'/bigwig_out.lzma'
        compress_pickle.dump(bigwig_out, output_filepath, compression="lzma")
        print('done, saving bigwig_out at: '+output_filepath)
    else:
        print('done, returning bigwig_out as object')
        return bigwig_out