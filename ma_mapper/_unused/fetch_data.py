#%%
import sys
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
def extract_bam(bam_file, chrom, start_list, end_list, strand, offset = 5, probe_length =100, smoothing_length= 100):
    #print(chrom, start_list, end_list, strand)
    import pysam
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
    import pysam
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
        logging.info('done, saving bam_min at: '+output_filepath_min)
        output_filepath_max = output_dir+'/bam_max.lzma'
        compress_pickle.dump(bam_max, output_filepath_max, compression="lzma")
        logging.info('done, saving bam_max at: '+output_filepath_max)
        output_filepath_forward = output_dir+'/bam_forward.lzma'
        compress_pickle.dump(bam_forward, output_filepath_forward, compression="lzma")
        logging.info('done, saving bam_forward at: '+output_filepath_forward)
        output_filepath_reverse = output_dir+'/bam_reverse.lzma'
        compress_pickle.dump(bam_reverse, output_filepath_reverse, compression="lzma")
        logging.info('done, saving bam_reverse at: '+output_filepath_reverse)
    else:
        logging.info('done, returning bam_min, bam_max, bam_forward, bam_reverse as objects')
        return bam_min, bam_max, bam_forward, bam_reverse
#%%
def extract_vcf(vcf_file, chrom, start_list, end_list, strand, query_key):
    windows = []
    for i in range(len(start_list)):
        # zero-based to one-based conversion [BED -> VCF]
        start = start_list[i] +1
        end = end_list[i]
        genome_position = chrom + ":" + str(start) + "-" + str(end)
        window = np.zeros(end-start +1, dtype=object) #offset the +1 at the start
        import cyvcf2
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
    logging.info('extract from vcf target: '+ vcf_input)
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
    import cyvcf2
    with ProcessPoolExecutor(max_workers=40) as executor:
        results  = executor.map(extract_vcf, vcf_call_list, chrom_list, start_list,end_list,strand_list, repeat(query_key))
    vcf_out = []
    for result in results:
        vcf_out.append(result)
    if save_to_file == True:
        import compress_pickle
        output_filepath = output_dir+'/vcf_out.lzma'
        compress_pickle.dump(vcf_out, output_filepath, compression="lzma")
        logging.info('done, saving vcf_out at: '+output_filepath)
    else:
        logging.info('done, returning vcf_out as object')
        return vcf_out
#%%
if sys.version_info >= (3, 8, 0):
    from typing import Literal
else:
    from typing_extensions import Literal
_AGEARG = Literal[None,'extract','calibrate']
_COUNTARG = Literal['human_ref','coverage','common', 'common_raw', 'total_raw']
def extract_maf(maf_file, chrom, start, end, strand,target_species = "Homo_sapiens", count_arg:_COUNTARG = 'human_ref',age=None, age_arg:_AGEARG = None, age_table_file = None):
    #print(target_species,chrom, start, end, strand)
    if age_arg is not None:
        if isinstance(age_table_file, str):
            if (os.path.isfile(age_table_file) == True):
                species = pd.read_table(age_table_file, sep='\t')
            else:
                print('metadata file not found')
        else:
            species = age_table_file
    maf_id = target_species+'.'+chrom
    from Bio.AlignIO import MafIO
    index_maf = MafIO.MafIndex(maf_file+".mafindex", maf_file, maf_id) 
    if strand =='-':
        n_strand = -1
        results =index_maf.get_spliced(start,end,n_strand)
    else:
        n_strand = 1
        results =index_maf.get_spliced(start,end,n_strand)
    collector = {}
    for seqrec in results:
        if seqrec.id.split('.')[0] == target_species:  
            collector[seqrec.id.split('.')[0]]= np.array(list(seqrec.seq.lower()))
            age_measure = 0
            for seqrec in results:
                if seqrec.id.split('.')[0] != target_species:
                    if age_arg is None:
                        collector[seqrec.id.split('.')[0]]= np.array(list(seqrec.seq.lower()))
                    elif age_arg == 'calibrate':
                        species_age = species[species.meta_name == seqrec.id.split('.')[0]]['Estimated Time (MYA)'].to_list()[0]
                        if age >= species_age:
                            collector[seqrec.id.split('.')[0]]= np.array(list(seqrec.seq.lower()))
                    elif age_arg == 'extract':
                        species_age = species[species.meta_name == seqrec.id.split('.')[0]]['Estimated Time (MYA)'].to_list()[0]
                        if species_age >= age_measure:
                            age_measure = species_age
    if age_arg == 'extract':
        return age_measure
    else:
        array_transposed=np.array(list(collector.values())).transpose()
        output_array=[]
        for ref_pos in array_transposed:
            (unique, counts) = np.unique(ref_pos, return_counts=True)
            frequencies = dict(zip(unique, counts))
            ref_allele=ref_pos[0]
            #frequencies.pop('-', None) #->count deletion/insertion
            total = sum(frequencies.values())
            if count_arg == 'human_ref':
                alt_count = total - frequencies[ref_allele]
                alt_freq=alt_count/total
                output_array.append(alt_freq)
            elif count_arg == 'coverage':
                output_array.append(total)
            elif count_arg == 'common':
                common_allele = max(frequencies, key=frequencies.get)
                common_count = frequencies[common_allele]
                common_freq = common_count/total
                output_array.append(common_freq)
            elif count_arg == 'common_raw':
                common_allele = max(frequencies, key=frequencies.get)
                common_count = frequencies[common_allele]
                output_array.append(common_count)
            elif count_arg == 'total_raw':
                output_array.append(total)
        return output_array
#%%
def fetch_maf(metadata_input, maf_input,output_dir = None, separated_maf = False, target_species = 'Homo_sapiens', count_arg = 'human_ref', save_to_file = False, custom_id = False, age_arg:_AGEARG = None, age_table_file = None, age_list=None):
    if isinstance(metadata_input, str):
        if os.path.isfile(metadata_input):
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
    log_path = output_dir+'/'+'maf_extract.log'
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
    logging.info('extract from maf target: '+ maf_input)
    if custom_id == False:
        meta_id = 'sample_n'+ metadata.index.astype(str)
        metadata['meta_id'] = meta_id
    else:
        metadata['meta_id'] = metadata.iloc[:,4]
        meta_id = metadata.iloc[:,4].unique()
    maf_call_list = []
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
        if separated_maf == True:
            maf_file = maf_input + '.' + chrom
        else:
            maf_file = maf_input
        maf_call_list.append(maf_file)
    if age_arg is None:
        logging.info('normal maf extraction without calibration')
        with ProcessPoolExecutor(max_workers=40) as executor:
            results  = executor.map(extract_maf, maf_call_list, chrom_list, start_list, end_list, strand_list, repeat(target_species) ,repeat(count_arg))
    elif age_arg == 'extract':
        logging.info('extract age depth for further use')
        with ProcessPoolExecutor(max_workers=40) as executor:
            results  = executor.map(extract_maf, maf_call_list, chrom_list, start_list, end_list, strand_list, repeat(target_species) ,repeat(count_arg), repeat(age_list), repeat(age_arg), repeat(age_table_file))
    elif age_arg == 'calibrate':
        if age_list is None:
            logging.info('calibrate age: no age_list provided, fetching from metadata')
            age_list = []
            for uniq_meta_id in meta_id:
                age_list.append(metadata_by_id.te_age.unique()[0])
        else:
            logging.info('calibrate age: using preexisting age_list')
            age_list = age_list
        with ProcessPoolExecutor(max_workers=40) as executor:
            results  = executor.map(extract_maf, maf_call_list, chrom_list, start_list, end_list, strand_list, repeat(target_species) ,repeat(count_arg), age_list, repeat(age_arg), repeat(age_table_file))
    maf_out = []
    for result in results:
        maf_out.append(result)
    if save_to_file == True:
        import compress_pickle
        output_filepath = output_dir+'/maf_out.lzma'
        compress_pickle.dump(maf_out, output_filepath, compression="lzma")
        logging.info('done, saving maf_out at: '+output_filepath)
    else:
        logging.info('done, returning maf_out as object')
        return maf_out
#%%
def extract_bed(bed_file, chrom, start_list, end_list, strand):
    drawings = []
    import pybedtools
    if isinstance(bed_file, str):
        if (os.path.isfile(bed_file) == True):
            bed_file = pybedtools.BedTool(bed_file)
        else:
            print('bed file not found')
    else:
        bed_file = bed_file
    for i in range(len(start_list)):
        start = start_list[i]
        end = end_list[i]
        prep_dict = {'chrom':[chrom], 'start':[start], 'end':[end],'name':['temp'],'score':[10], strand:[strand]}
        temp_df = pd.DataFrame(prep_dict)
        metadata_bed=pybedtools.BedTool.from_dataframe(temp_df)
        temp_intersect=metadata_bed.intersect(bed_file,wa=True, wb =True, s=True)
        fragment_length = end - start
        canvas=np.zeros(fragment_length)
        for intersect in temp_intersect:
            start_intersect = int(intersect[7]) #thickEnd
            end_intersect = int(intersect[8]) #itemRgb
            score_intersect = float(intersect[10])
            if strand == '+':
                start_intersect_position = start_intersect - start
                end_intersect_postion = end_intersect - start
                #canvas[start_intersect_position:end_intersect_postion] += np.full(end_intersect-start_intersect, score_intersect)
            else:
                start_intersect_position = end - end_intersect +1
                end_intersect_postion = end - start_intersect +1
            canvas[start_intersect_position:end_intersect_postion] += score_intersect
        #if strand == '-':
            #canvas = np.flip(canvas)
        drawings.append(canvas)
    if strand == '-':
        drawings.reverse()
    bed_out = np.concatenate(drawings)
    return bed_out
#%%
def fetch_bed(metadata_input, bed_input, output_dir = None, save_to_file = False, custom_id = False):
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
    if isinstance(bed_input, str):
        if (os.path.isfile(bed_input) == True):
            import pybedtools
            bed_file = pybedtools.BedTool(bed_input)
        else:
            print('bed file not found')
    else:
        bed_file = bed_input
    import logging
    log_path = output_dir+'bed_extract.log'
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
    logging.info('extract from bam file: '+ bed_input)
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
        metadata_by_id = metadata[metadata.id == uniq_meta_id]
        chrom_list.append(metadata_by_id.iloc[:,0].unique()[0])
        start_list.append(metadata_by_id.iloc[:,1].to_list())
        end_list.append(metadata_by_id.iloc[:,2].to_list())
        strand_list.append(metadata_by_id.iloc[:,3].unique()[0])
    with ProcessPoolExecutor(max_workers=40) as executor:
        results  = executor.map(extract_bed, repeat(bed_file), chrom_list, start_list,end_list,strand_list)
    bed_out = []
    for result in results:
        bed_out.append(result)
    if save_to_file == True:
        import compress_pickle
        output_filepath = output_dir+'/bed_out.lzma'
        compress_pickle.dump(bed_out, output_filepath, compression="lzma")
        logging.info('done, saving bed_out at: '+output_filepath)
    else:
        logging.info('done, returning bed_out as object')
        return bed_out
#%%
def extract_bigwig(bigwig_file, chrom, start_list, end_list, strand):
    bigwig_arrays = []
    for i in range(len(start_list)):
        # zero-based to one-based conversion [BED -> VCF]
        start = start_list[i]
        end = end_list[i]
        import pyBigWig
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
    log_path = output_dir+'bigwig_extract.log'
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
    logging.info('extract from bigwig target: '+ bigwig_input)
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
        logging.info('done, saving bigwig_out at: '+output_filepath)
    else:
        logging.info('done, returning bigwig_out as object')
        return bigwig_out
#%%
def main():
    print('main process')

if __name__ == '__main__':
    main()
