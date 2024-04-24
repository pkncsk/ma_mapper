#%%
import math
import numpy as np
import pandas as pd
import pysam
import os
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
def extract_bam(bam_file, chrom, start_list, end_list, strand, offset = 5, probe_length =100, smoothing_length= 100):
    #print(chrom, start_list, end_list, strand)
    if isinstance(bam_file, str):
        bam_file = pysam.AlignmentFile(bam_file, "rb")
    normal = normal_array(width=probe_length, sigma=smoothing_length, odd=1)
    result_min = []
    result_max = []
    result_forward = []
    result_reverse = []
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

        output_min = np.minimum(profile_normals_forward, profile_normals_reverse)
        output_max = np.maximum(profile_normals_forward, profile_normals_reverse)
        output_forward = np.array(profile_reads_forward)
        output_reverse = np.array(profile_reads_reverse)

        if strand == '-':
            output_min = np.flip(output_min)
            output_max = np.flip(output_max)
            output_forward = np.flip(output_forward)
            output_reverse = np.flip(output_reverse)

        output_min = output_min[probe_length:probe_length+fragment_length]
        output_max = output_max[probe_length:probe_length+fragment_length]
        output_forward = output_forward[probe_length:probe_length+fragment_length]
        output_reverse = output_reverse[probe_length:probe_length+fragment_length]
    
        result_min.append(output_min)
        result_max.append(output_max)
        result_forward.append(output_forward)
        result_reverse.append(output_reverse)
    #failsafe in case of no mapped reads
    if len(result_min) != 0:
        if strand == '-':
            result_min.reverse()
        np_min = np.concatenate(result_min)
    else:
        np_min = np.zeros(np.sum(np.subtract(end_list,start_list)), dtype=np.uint8)
    
    if len(result_max) != 0:
        if strand == '-':
            result_max.reverse()
        np_max = np.concatenate(result_max)
    else:
        np_max = np.zeros(np.sum(np.subtract(end_list,start_list)), dtype=np.uint8)
    
    if len(result_forward) != 0:
        if strand == '-':
            result_forward.reverse()
        np_forward = np.concatenate(result_forward)
    else:
        np_forward = np.zeros(np.sum(np.subtract(end_list,start_list)), dtype=np.uint8)

    if len(result_reverse) != 0:
        if strand == '-':
            result_reverse.reverse()
        np_reverse = np.concatenate(result_reverse)
    else:
        np_reverse = np.zeros(np.sum(np.subtract(end_list,start_list)), dtype=np.uint8)
    return np_min, np_max, np_forward, np_reverse
#%%
def fetch_bam(metadata_input, bam_input, output_dir = None, offset = 5, probe_length =100, smoothing_length= 100, save_to_file = False, custom_id = False):
    bam_file = pysam.AlignmentFile(bam_input, "rb")
    global normal
    normal = normal_array(width=probe_length, sigma=smoothing_length, odd=1)
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
    if custom_id == False:
        meta_id = 'sample_n'+ metadata.index.astype(str)
        metadata['meta_id'] = meta_id
    else:
        metadata['meta_id'] = metadata.iloc[:,4]
        meta_id = metadata.iloc[:,4].unique()   
    bam_min = list()
    bam_max = list()
    bam_forward = list()
    bam_reverse = list()
    for uniq_meta_id in meta_id:
        metadata_by_id = metadata[metadata.meta_id == uniq_meta_id]
        chrom = metadata_by_id.iloc[:,0].unique()[0]
        start_list = metadata_by_id.iloc[:,1].to_list()
        end_list = metadata_by_id.iloc[:,2].to_list()
        strand = metadata_by_id.iloc[:,3].unique()[0]
        np_min, np_max, np_forward, np_reverse = extract_bam(bam_file,chrom, start_list, end_list, strand, offset, probe_length, smoothing_length)
        bam_min.append(np_min)
        bam_max.append(np_max)
        bam_forward.append(np_forward)
        bam_reverse.append(np_reverse)
    import compress_pickle
    if save_to_file == True:
        output_filepath_min = output_dir+'/bam_min.lzma'
        compress_pickle.dump(bam_min, output_filepath_min, compression="lzma")
        print('done, saving bam_min at: '+output_filepath_min)
        output_filepath_max = output_dir+'/bam_max.lzma'
        compress_pickle.dump(bam_max, output_filepath_max, compression="lzma")
        print('done, saving bam_max at: '+output_filepath_max)
        output_filepath_forward = output_dir+'/bam_forward.lzma'
        compress_pickle.dump(bam_forward, output_filepath_forward, compression="lzma")
        print('done, saving bam_forward at: '+output_filepath_forward)
        output_filepath_reverse = output_dir+'/bam_reverse.lzma'
        compress_pickle.dump(bam_reverse, output_filepath_reverse, compression="lzma")
        print('done, saving bam_reverse at: '+output_filepath_reverse)
    else:
        print('done, returning bam_min, bam_max, bam_forward, bam_reverse as objects')
        return bam_min, bam_max, bam_forward, bam_reverse