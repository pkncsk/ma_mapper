import numpy as np
import pandas as pd
import cyvcf2
import os
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
#%%
def extract_vcf(vcf_file, chrom, start_list, end_list, strand, query_key):
    windows = []
    for i in range(len(start_list)):
        # zero-based to one-based conversion [BED -> VCF]
        start = start_list[i] +1
        end = end_list[i]
        genome_position = chrom + ":" + str(start) + "-" + str(end)
        window = np.zeros(end-start +1, dtype=object) #offset the +1 at the start
        try:
            vcf = cyvcf2.VCF(vcf_file)
        except OSError:
                print("cannot find "+genome_position)
        for variant in vcf(genome_position):
            if variant.var_type == 'snp':            
                relative_position = variant.POS-start
                try:
                    window[relative_position] +=  variant.INFO[query_key]
                except KeyError:
                    window[relative_position] +=  0
            if strand == '-':
                window = np.flip(window)
        windows.append(window)
    if strand == '-':
        windows.reverse()
    window_out = np.concatenate(windows)
    return window_out
#%%
def fetch_vcf(metadata_input, vcf_input, query_key = 'AF', output_dir = None, vcf_format = None, save_to_file = False, custom_id = False):
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
    log_path = output_dir+'vcf_extract.log'
    #setup logger
    logging.root.handlers = []
    logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    handlers = [
                        logging.FileHandler(log_path, mode = 'a'),
                        logging.StreamHandler()
                        ]
                    )
    print('extract from vcf target: '+ vcf_input)
    if custom_id == False:
        meta_id = 'sample_n'+ metadata.index.astype(str)
        metadata['meta_id'] = meta_id
    else:
        metadata['meta_id'] = metadata.iloc[:,4]
        meta_id = metadata.iloc[:,4].unique()    
    vcf_call_list = []
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
        if vcf_format == 'gnomad':
            vcf_file = vcf_input + 'gnomad.genomes.v3.1.1.sites.'+chrom+'.vcf.gz'
        else:
            vcf_file = vcf_input
        vcf_call_list.append(vcf_file)
    with ProcessPoolExecutor(max_workers=40) as executor:
        results  = executor.map(extract_vcf, vcf_call_list, chrom_list, start_list,end_list,strand_list, repeat(query_key))
    vcf_out = []
    for result in results:
        vcf_out.append(result)
    if save_to_file == True:
        import compress_pickle
        output_filepath = output_dir+'/vcf_out.lzma'
        compress_pickle.dump(vcf_out, output_filepath, compression="lzma")
        print('done, saving vcf_out at: '+output_filepath)
    else:
        print('done, returning vcf_out as object')
        return vcf_out
#%%