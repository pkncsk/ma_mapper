#%%
import pysam
import math
import numpy as np
import pandas as pd
import os
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
#%%
#-> BAM file overlay 
def normal_array(width=1, sigma=1, odd=0):
    ''' Returns an array of the normal distribution of the specified width '''
    sigma2 = float(sigma) ** 2
    def normal_func(x):
        return math.exp(-x * x / (2 * sigma2))
    # width is the half of the distribution
    values = list(map(normal_func, range(-width, width + odd)))
    values = np.array(values)
    return values
#%%
def extract_bam(bam_file, chrom, start, end, strand, offset = 5, probe_length =100, smoothing_length= 100):

    data_min = list()
    data_max = list()
    data_forward = list()
    data_reverse = list()
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

        #if mode == 'smooth-min':
        output_min = np.minimum(profile_normals_forward, profile_normals_reverse)
        #elif mode == 'smooth-max':
        output_max = np.maximum(profile_normals_forward, profile_normals_reverse)
        #elif mode == 'forward-read':
        output_forward = np.array(profile_reads_forward)
        #elif mode == 'reverse-read':
        output_reverse = np.array(profile_reads_reverse)
        #else:
        #    output_signal = np.maximum(profile_reads_forward, profile_reads_reverse)

        if strand == '-':
            output_min = np.flip(output_min)
            output_max = np.flip(output_max)
            output_forward = np.flip(output_forward)
            output_reverse = np.flip(output_reverse)

        output_min = output_min[probe_length:probe_length+fragment_length]
        output_max = output_max[probe_length:probe_length+fragment_length]
        output_forward = output_forward[probe_length:probe_length+fragment_length]
        output_reverse = output_reverse[probe_length:probe_length+fragment_length]
     
        
        data_min.extend(output_min)
        data_max.extend(output_max)
        data_forward.extend(output_forward)
        data_reverse.extend(output_reverse)
        
    #print(data)
    np_min = np.array(data_min)
    np_max = np.array(data_max)
    np_forward = np.array(data_forward)
    np_reverse = np.array(data_reverse)
    
    return np_min, np_max, np_forward, np_reverse
#%%
def fetch_bam(metadata_input, bam_input, output_dir = None, offset = 5, probe_length =100, smoothing_length= 100):
    bam_file = pysam.AlignmentFile(bam_input, "rb")
    global normal
    normal = normal_array(width=probe_length, sigma=smoothing_length, odd=1)
    if (os.path.isfile(metadata_input) == True):
        metadata = pd.read_csv(metadata_input, delim_whitespace=True)
    else:
        metadata = metadata_input
    if output_dir is None:
        output_dir = '/'.join(str.split(metadata_input, sep ='/')[:-1])
    import logging
    log_path = output_dir+'bam_extract.log'
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
    logging.info('extract from bam file: '+ bam_input)
    bam_min = list()
    bam_max = list()
    bam_forward = list()
    bam_reverse = list()
    for idx, row in metadata.iterrows():
        genoname = row.genoName
        genostart = row.genoStart
        genoend = row.genoEnd
        strand = row.strand
        np_min, np_max, np_forward, np_reverse = extract_bam(bam_file ,genoname, genostart, genoend, strand, offset, probe_length, smoothing_length)
        bam_min.append(np_min)
        bam_max.append(np_max)
        bam_forward.append(np_forward)
        bam_reverse.append(np_reverse)
    import compress_pickle
    output_filepath_min = output_dir+'/bam_min.lzma'
    compress_pickle.dump(bam_min, output_filepath_min, compression="lzma")
    logging.info('done, saving te bam_min at: '+output_filepath_min)
    output_filepath_max = output_dir+'/bam_max.lzma'
    compress_pickle.dump(bam_max, output_filepath_max, compression="lzma")
    logging.info('done, saving te bam_max at: '+output_filepath_max)
    output_filepath_forward = output_dir+'/bam_forward.lzma'
    compress_pickle.dump(bam_forward, output_filepath_forward, compression="lzma")
    logging.info('done, saving te bam_forward at: '+output_filepath_forward)
    output_filepath_reverse = output_dir+'/bam_reverse.lzma'
    compress_pickle.dump(bam_reverse, output_filepath_reverse, compression="lzma")
    logging.info('done, saving te bam_reverse at: '+output_filepath_reverse)
#%%
def main():
    print('main process')

if __name__ == '__main__':
    main()