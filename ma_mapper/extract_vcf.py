import numpy as np
import pandas as pd
import cyvcf2
import os
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from . import logger
#%%
def extract_vcf(vcf_file, chrom, start_list, end_list, strand, query_key):
    windows = []
    for i in range(len(start_list)):
        # zero-based to one-based conversion [BED -> VCF]
        start = start_list[i] +1
        end = end_list[i]
        genome_position = f'{chrom}:{start}-{end}'
        window = np.zeros(end-start +1, dtype=float) #offset the +1 at the start
        try:
            vcf = cyvcf2.VCF(vcf_file)
        except OSError:
                logger.error("cannot find ",genome_position)
        for variant in vcf(genome_position):
            if variant.var_type == 'snp':            
                relative_position = variant.POS-start
                try:
                    window[relative_position] +=  variant.INFO[query_key]
                except KeyError:
                    window[relative_position] +=  0.0
            if strand == '-':
                window = np.flip(window)
        windows.append(window)
    if strand == '-':
        windows.reverse()
    window_out = np.concatenate(windows)
    return window_out.astype(float)
#%%
def vcf_io(metadata, vcf, query_key = 'AF', output_dir = None, vcf_format = None, save_to_file = False, custom_id = False, custom_prefix='entry'):
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

    grouped = metadata_local.groupby('meta_id', sort=False)
    chrom_list = grouped.apply(lambda x: x.iloc[:,0].unique()[0], include_groups=False).tolist()
    start_list = grouped.apply(lambda x: x.iloc[:,1].tolist(), include_groups=False).tolist()
    end_list = grouped.apply(lambda x: x.iloc[:,2].tolist(), include_groups=False).tolist()
    strand_list = grouped.apply(lambda x: x.iloc[:,5].unique()[0], include_groups=False).tolist()
    vcf_call_list = []
    for chrom in chrom_list:
        if vcf_format == 'gnomad':
            vcf_file = f'{vcf}/gnomad.genomes.v3.1.1.sites.{chrom}.vcf.gz'
        else:
            vcf_file = vcf
        vcf_call_list.append(vcf_file)
    with ProcessPoolExecutor(max_workers=40) as executor:
        results  = executor.map(extract_vcf, vcf_call_list, chrom_list, start_list,end_list,strand_list, repeat(query_key))
    vcf_out = []
    for result in results:
        vcf_out.append(result)
    if save_to_file == True:
        import compress_pickle
        output_filepath = f'{output_dir}/vcf_out.lzma'
        compress_pickle.dump(vcf_out, output_filepath, compression="lzma")
        logger.info('done, saving vcf_out at: ', output_filepath)
    else:
        logger.info('done, returning vcf_out as object')
        return vcf_out
#%%