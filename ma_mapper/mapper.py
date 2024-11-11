#%%
from turtle import back
import pandas as pd

from ma_mapper import extract_vcf
from . import sequence_alignment
import numpy as np
import sys
import os
from . import logger
import collections
#from typing import Literal <- python3.8+
if sys.version_info >= (3, 8, 0):
    from typing import Literal, Tuple, List
else:
    from typing_extensions import Literal, Tuple, List
#%%
def load_alignment_file(alignment_file):
    if isinstance(alignment_file, str) == True:
        from Bio import SeqIO
        records = (r for r in SeqIO.parse(alignment_file, "fasta"))
    else:
        records = alignment_file
    return records

def extract_metadata_from_alignment(alignment_file, save_to_file = False, output_file =None):
    records = load_alignment_file(alignment_file)
    meta_id_list = []; chrom_list = []; start_list = []; end_list = []; strand_list = []; score_list= []
    import re
    for record in records:
        identifier=record.name
        uniq_meta_id, chrom, start, end, strand = re.match(r'([^:]+)::([^:]+):(\d+)-(\d+)\((.)\)', identifier).groups()
        chrom_list.append(chrom)
        start_list.append(int(start))
        end_list.append(int(end))
        meta_id_list.append(uniq_meta_id)
        score_list.append(20)
        strand_list.append(strand)
    new_metadata_dict = {'chrom':chrom_list,'start':start_list,'end':end_list,'name':meta_id_list,'score':score_list,'strand':strand_list}
    new_metadata = pd.DataFrame(new_metadata_dict)
    if save_to_file == True:
        if output_file is not None:
            output_filepath = output_file
        else:
            import os
            output_filepath =f'{os.path.dirname(os.path.abspath(__file__))}/metadata_aligned.txt'
        new_metadata.to_csv(output_filepath, sep='\t', index= False)
        logger.info('save new metadata at', output_filepath)
    else:
        return new_metadata

def parse_alignment(alignment_file:str, 
                    save_to_file: bool =False, 
                    output_file: str| None=None):
    if output_file is None:
        if isinstance(alignment_file, str):
            output_file = f'{alignment_file}.parsed'
        else:
            output_file = f'{os.path.dirname(os.path.abspath(__file__))}/alignment.parsed'
    output_dir = '/'.join(str.split(output_file, sep ='/')[:-1])
    logger.info('parse alignment')
    if isinstance(alignment_file, str):
        if (os.path.isfile(alignment_file) == True):
            records = load_alignment_file(alignment_file)
        else:
            logger.error(f"Alignment file '{alignment_file}' not found.")
            raise FileNotFoundError(f"Alignment file '{alignment_file}' not found.")
    else:
        records = alignment_file
    #counting sequence records
    seq_count = 0
    for i, value in enumerate(records):
        if seq_count == 0:
            seq_length = len(value)
        seq_count = seq_count + 1
    if seq_count == 0:
        logger.info('no records')
    else:
        seq_id_list = list()
        for idx, content in enumerate(records):
            seq_id_list.append(content.name)
    #create parsed array
    parsed_array = np.zeros((seq_count, seq_length), dtype=np.uint8)
    if isinstance(alignment_file, str):
        if (os.path.isfile(alignment_file) == True):
            records = load_alignment_file(alignment_file)
        else:
            logger.error('alignment file not found')
    else:
        records = alignment_file
    for i, value in enumerate(records):
        parsed_array[i, :] = np.array(str(value.seq).upper(), 'c').view(np.uint8)
    
    hash_table = {45:0, 65:1, 67:2, 84:3, 71:4, 78:5}
    result_array = np.ndarray(parsed_array.shape)
    for k in hash_table:
        result_array[parsed_array == k] = hash_table[k]
    if save_to_file == True:
        import h5py
        with h5py.File(output_file, 'w') as hf:
            hf.create_dataset("parsed_array",  data=result_array, compression="lzf")
        logger.info('done, saving parsed array at: ',output_file)
    else: 
        return result_array
#%%
def import_parsed_array(parsed_array_file):
    import h5py
    with h5py.File(f'{import_parsed_array}/parsed_array.h5', 'r') as hf:
        parsed_array = hf['result_array'][:]  
    return parsed_array
#%%#FIXME: find a better way to sort your alignment (and why does it need to be sorted?)
def sort_parsed_array(parsed_array, aligned_seqeunce):
    if isinstance(parsed_array, str) == True:
        parsed_array = import_parsed_array(parsed_array)
    if isinstance(aligned_seqeunce, str) == True:
        from Bio import SeqIO
        records = (r for r in SeqIO.parse(aligned_seqeunce, "fasta"))
    else:
        records = aligned_seqeunce
    align_order = []
    seq_id_list = []
    seq_row_nmber_list=[]
    for i, value in enumerate(records):
        #print(i)
        seq_id, seq_row_nmber=value.name.split(sep='(')[0].split(sep='.')
        seq_id_list.append(int(seq_id))
        seq_row_nmber_list.append(int(seq_row_nmber))
        align_order.append(i)

    seq_id_list_sorted, seq_row_nmber_list_sorted, parsed_array_sorted,align_order_sort = zip(*sorted(zip(seq_id_list, seq_row_nmber_list, parsed_array,align_order)))
    parsed_array_sorted = np.array(parsed_array_sorted)
    return parsed_array_sorted
#%%
def create_filter(parsed_array, col_threshold = 0.50, col_content_threshold = 0.10, row_threshold = 0.50, save_to_file = False, output_file = None):
    if isinstance(parsed_array, str) == True:
        parsed_array = import_parsed_array(parsed_array)
    col_counts = np.zeros(parsed_array.shape[1])
    col_counts_max = np.zeros(parsed_array.shape[1])  
    for i, row in enumerate(parsed_array):
        nonzeros = np.nonzero(row)[0]
        col_counts[nonzeros] += 1
        col_counts_max[nonzeros[0]:nonzeros[-1]] += 1
    col_filter = (col_counts > col_counts_max * col_threshold) & (col_counts >= parsed_array.shape[0] * col_content_threshold)
    row_ratio = np.zeros(parsed_array.shape[0])
    for i, row in enumerate(parsed_array[:, col_filter]):
        nonzeros = np.nonzero(row)[0]
        if nonzeros.size and float(nonzeros[-1] - nonzeros[0]) > 0:
            row_ratio[i] = len(nonzeros) / (float(nonzeros[-1] - nonzeros[0]))
            #print(row_ratio[i])
    row_filter = row_ratio >= row_threshold
    for i, row in  enumerate(parsed_array[row_filter, :]):
        nonzeros = np.nonzero(row)[0]

        col_counts[nonzeros] += 1
        col_counts_max[nonzeros[0]:nonzeros[-1]] += 1
    col_filter = (col_counts > col_counts_max * col_threshold) & (col_counts >= parsed_array.shape[0] * col_content_threshold)
    filters = [row_filter,col_filter]
    if save_to_file == True:
        import pickle
        if output_file is None:
            output_file = 'filter.p'
        pickle.dump(filters, open(output_file, "wb" ))
    else:
        return filters
#%%
def create_filter_failed(parsed_array, col_threshold = 0.50, row_threshold = 0.50, save_to_file = False, output_file = None):
    if isinstance(parsed_array, str) == True:
        parsed_array = import_parsed_array(parsed_array)
    col_filter = list(range(parsed_array.shape[1]))
    row_filter = list(range(parsed_array.shape[0]))
    while True:
        # Filter columns
        new_col_filter = []
        for idx in col_filter:
            nonzero_ratio = np.count_nonzero(parsed_array[row_filter, idx]) / len(row_filter)
            if nonzero_ratio > col_threshold:
                new_col_filter.append(idx)
        # Filter rows
        new_row_filter = []
        for idx in row_filter:
            nonzero_ratio = np.count_nonzero(parsed_array[idx, new_col_filter]) / len(new_col_filter)
            if nonzero_ratio > row_threshold:
                new_row_filter.append(idx)
        
        # Check for convergence
        if new_col_filter == col_filter and new_row_filter == row_filter:
            break
        col_filter, row_filter = new_col_filter, new_row_filter    

    filters = [row_filter,col_filter]
    if save_to_file == True:
        import pickle
        if output_file is None:
            output_file = 'filter.p'
        pickle.dump(filters, open(output_file, "wb" ))
    else:
        return filters
#%%
def map_data(data_file:str,
             sorted_parsed_array: np.ndarray, 
             filters: bool = True,
             custom_filters = None,
             col_threshold: int = 0.50, 
             col_content_threshold:int = 0.10, 
             row_threshold:int = 0.50,
             nested_data:bool=None)-> np.ndarray:
    if isinstance(data_file, str) == True:
        import compress_pickle
        data_file = compress_pickle.load(data_file)
    def get_unique_element_types(matrix):
        unique_types = set()
        for row in matrix:
            for element in row:
                unique_types.add(type(element))
        return unique_types
    
    if nested_data is None:
        unique_types = get_unique_element_types(data_file)
        logger.info(unique_types)
        if any(t in unique_types for t in (int,float,np.uint8,np.uint16, np.float32, np.float64)) :
            nested_data = False
        #elif any(t in unique_types for t in (str,set(),collections.Counter, list, dict, pd.DataFrame)):
        #    nested_data = True
        else:
            nested_data = True
    #def check_list_type(data, type_check):
    #    return all(isinstance(item, type_check) for sublist in data for item in sublist)
    #if check_list_type(sorted_parsed_array, (int, float, type(np.nan))):
    #    canvas = np.zeros(sorted_parsed_array.shape, dtype=object)
    #else:
    logger.info(f'nested_data:{nested_data}')
    if nested_data:
        canvas = np.full(sorted_parsed_array.shape, None)
    else:
        canvas = np.zeros(sorted_parsed_array.shape, dtype=float)
        
    for idx, row in enumerate(sorted_parsed_array):
        #print(f"data_file[{idx}]: {data_file[idx]}")
        #print(f"canvas[{idx}]: {canvas[idx, row>0]}")
        canvas[idx, row>0] = data_file[idx]
    mapped_data = canvas
    if filters:
        if custom_filters is None: 
            filters=create_filter(sorted_parsed_array, col_threshold = col_threshold, col_content_threshold = col_content_threshold, row_threshold = row_threshold)
        else:
            filters = custom_filters
        row_filter = filters[0]
        col_filter = filters[1]
        mapped_data=mapped_data[np.ix_(row_filter,col_filter)]
    return mapped_data
#%%
def base_count(alignment:pd.DataFrame):
    alignment_transposed = alignment.transpose()
    output_array = []
    for idx, ref_count in enumerate(alignment_transposed):
        (unique, counts) = np.unique(ref_count, return_counts=True)
        frequencies = dict(zip(unique, counts))
        base_count = []
        for base in [0,1,2,3,4]:
            nucl_count = frequencies.get(base,0)
            base_count.append(nucl_count)
        output_array.append(base_count)
    base_count=np.array(output_array).transpose()
    return base_count
#FIXED: wrong nonzero
_NORM_METHOD = Literal['average','perc_coverage','median']
def normalise(alignment: str|pd.DataFrame,
              mapped_data: pd.DataFrame, 
              method:_NORM_METHOD = 'average', 
              ):
    if method == 'average':
        normalizer = np.count_nonzero(alignment, axis=0)
        data_sum = np.nansum(mapped_data, axis=0)
        normalised_data = data_sum/normalizer 
    if method == 'perc_coverage':
        normalizer = np.count_nonzero(alignment, axis=0)
        data_count = np.count_nonzero(mapped_data, axis=0)
        normalised_data = data_count/normalizer * 100
    if method == 'median':
        normalised_data = np.nanmedian(mapped_data, axis = 0)
    return normalised_data

def flank_sequence_io(metadata: pd.DataFrame,
                    source_fasta: str,
                    metadata_out:bool = False, 
                    extension_length:int = 500) -> Tuple[pd.DataFrame, pd.DataFrame]:

    front_metadata = metadata.reset_index(drop=True).copy()
    front_metadata.end = front_metadata.start.astype(int)
    front_metadata.start = front_metadata.start.astype(int) - extension_length
    
    back_metadata = metadata.reset_index(drop=True).copy()
    back_metadata.start = back_metadata.end.astype(int) 
    back_metadata.end = back_metadata.end.astype(int) + extension_length
    if metadata_out:
        return front_metadata, back_metadata
    else:
        front_seq = sequence_alignment.sequence_io(front_metadata, source_fasta, custom_id=False)
        back_seq = sequence_alignment.sequence_io(front_metadata, source_fasta, custom_id=False)
        
        front_list = []
        back_list = []

        for i, row in metadata.reset_index().iterrows():
            if row.strand == '+':
                front_list.append(front_seq[i])
                back_list.append(back_seq[i])
            else:
                front_list.append(back_seq[i])
                back_list.append(front_seq[i])
        front_parsed=parse_alignment(front_list)
        back_parsed =parse_alignment(back_list)
        return front_parsed, back_parsed
#steamline functions
def parse_and_filter(alignment_file: str,
                     filters:bool = True,
                     extension_length:int | None = None,
                     preprocess_out:bool = False,
                     **kwargs)-> Tuple[np.ndarray, pd.DataFrame]:
    
    aligned_parsed = parse_alignment(alignment_file)
    metadata_aligned = extract_metadata_from_alignment(alignment_file)
    metadata_aligned['original_order'] = metadata_aligned.index
    if filters:
        cf_kwargs = {key: value for key, value in kwargs.items() if key in ('col_threshold', 'row_threshold')}
        filters=create_filter(aligned_parsed, **cf_kwargs)
        row_filter = filters[0]
        col_filter = filters[1]
        aligned_filtered=aligned_parsed[np.ix_(row_filter,col_filter)]
        metadata_aligned_filtered=metadata_aligned.iloc[row_filter,:]
    else:
        aligned_filtered = aligned_parsed
        metadata_aligned_filtered = metadata_aligned
    #print(metadata_aligned_filtered.dtypes)
    if extension_length is not None:
        ffs_kwargs = {key: value for key, value in kwargs.items() if key in ('source_fasta')}
        front_parsed, back_parsed= flank_sequence_io(metadata = metadata_aligned_filtered, extension_length=extension_length, **ffs_kwargs)
        alignment_output = np.hstack((front_parsed, aligned_filtered, back_parsed))
    else:
        alignment_output = aligned_filtered
    
    if preprocess_out:
        return aligned_parsed, metadata_aligned, filters
    else:
        return alignment_output, metadata_aligned_filtered

def match_age_to_id_metadata(metadata: pd.DataFrame|str, 
                             age_table:str|None = None) -> pd.DataFrame:
    if isinstance(metadata, str):
        if os.path.isfile(metadata):
            metadata_df = pd.read_csv(metadata, sep='\t')
        else:
            logger.error('metadata file not found')
    elif isinstance(metadata, pd.DataFrame):
        metadata_df = metadata
    if isinstance(age_table,list):
        age_df_list = []
        for age_tbl in age_table:
            age_df_list.append(pd.read_csv(age_tbl, sep='\t'))
        age_df=pd.concat(age_df_list)
    elif isinstance(age_table, str):
        if os.path.isfile(age_table):
            age_df = pd.read_csv(age_table, sep='\t')
        else:
            logger.error('metadata file not found')
            raise FileNotFoundError(f"metadata file not found")
    elif isinstance(age_table, pd.DataFrame):
        age_df = age_table
    #assume BED format:
    bed_headers = ['chrom','start','end','name','score','strand']
    if metadata_df.shape[1] == 6:
        metadata_df.columns = bed_headers
    elif metadata_df.shape[1] > 6:
        current_columns = metadata_df.columns.tolist()
        metadata_df.columns = bed_headers + current_columns[6:]
    merged_table=metadata_df.merge(age_df, left_on='name', right_on='internal_id', how='left')
    merged_table=merged_table.drop(columns=['internal_id'])

    return merged_table

_FORMAT = Literal['bam_max','bam_min','bam_forward','bam_reverse','bigwig','bed']
def map_and_overlay(aligment:str,
            metadata:str,
            data_file:str,
            data_format:_FORMAT,
            filters:bool = True,
            extension_length: int|None = None,
            source_fasta:str = None,
            custom_id:bool=False,
            *args,**kwargs,
            )->Tuple[np.ndarray,pd.DataFrame,np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    fs_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('fs_')}
    pf_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('pf_')}
    md_kwargs = {key[3:]: value for key, value in kwargs.items() if key.startswith('md_')}
    kwargs = {key: value for key, value in kwargs.items() if not key.startswith(('pf_','md_','fs_'))}
    alignment_parsed, metadata_aligned, filters  = parse_and_filter(alignment_file=aligment, preprocess_out=True, **pf_kwargs)
    
    if data_format in ['read_max','read_min','read_forward','read_reverse', 'normal_max','normal_min','normal_forward','normal_reverse','read_sum']:
        from . import extract_bam
        output=extract_bam.bam_io(metadata= metadata, bam_file=data_file,bam_format=data_format, custom_id=custom_id,**kwargs)
    elif data_format in ['bigwig']:
        from . import extract_bigwig
        output=extract_bigwig.bigwig_io(metadata=metadata, bigwig=data_file, custom_id=custom_id,**kwargs)
    elif data_format in ['bed']:
        from . import extract_bed
        output=extract_bed.bed_io(metadata=metadata, bed=data_file, custom_id=custom_id,**kwargs)
    elif data_format in ['maf']:
        from . import extract_maf
        output=extract_maf.maf_io(metadata=metadata, maf=data_file,custom_id=custom_id,**kwargs)
    elif data_format in ['vcf']:
        output=extract_vcf.vcf_io(metadata=metadata, vcf=data_file, custom_id=custom_id, **kwargs)
    else:
        logger.error('dataformat field unspecified')
    if isinstance(metadata, str):
        if (os.path.isfile(metadata) == True):
            metadata_df = pd.read_csv(metadata, sep='\t', header=None)
        else:
            logger.error('metadata file not found')
    else:
        metadata_df = metadata   
    source_order = metadata_df.iloc[:,3].unique()
    output_sorted = []
    for _, row in metadata_aligned.iterrows():
        #print(row)
        #print(np.where(source_order == row['name']),row['name'])
        output_sorted.append(output[np.where(source_order == row['name'])[0][0]])
    
    overlay = map_data(data_file=output_sorted, sorted_parsed_array = alignment_parsed, custom_filters=filters, **md_kwargs)
    row_filter = filters[0]
    metadata_filtered=metadata_aligned.iloc[row_filter,:]

    if extension_length is not None:
        
        front_metadata, back_metadata = flank_sequence_io(metadata_filtered, metadata_out=True, extension_length=extension_length,source_fasta=source_fasta, **fs_kwargs)
        if data_format in ['read_max','read_min','read_forward','read_reverse', 'normal_max','normal_min','normal_forward','normal_reverse','read_sum']:
            front_output=extract_bam.bam_io(metadata= front_metadata, bam_file=data_file,bam_format=data_format, **kwargs)
            back_output=extract_bam.bam_io(metadata= back_metadata, bam_file=data_file,bam_format=data_format, **kwargs)
        elif data_format in ['bigwig']:
            front_output=extract_bigwig.bigwig_io(metadata=front_metadata, bigwig=data_file,custom_id=True,**kwargs)
            back_output=extract_bigwig.bigwig_io(metadata=back_metadata, bigwig=data_file,custom_id=True,**kwargs)
        elif data_format in ['bed']:
            front_output=extract_bed.bed_io(metadata=front_metadata, bed=data_file,custom_id=True,**kwargs)
            back_output=extract_bed.bed_io(metadata=back_metadata, bed=data_file,custom_id=True,**kwargs)
        elif data_format in ['maf']:
            front_output=extract_maf.maf_io(metadata=front_metadata, maf=data_file,custom_id=True,**kwargs)
            back_output=extract_maf.maf_io(metadata=back_metadata, maf=data_file,custom_id=True,**kwargs)
        elif data_format in ['vcf']:
            front_output=extract_vcf.vcf_io(metadata=front_metadata, vcf=data_file,custom_id=True,**kwargs)
            back_output=extract_vcf.vcf_io(metadata=back_metadata, vcf=data_file,custom_id=True,**kwargs)
        front_overlay = []
        back_overlay = []
        print(metadata_filtered.shape, len(front_output), len(back_output))
        for i, row in metadata_filtered.reset_index().iterrows():
            #print(i)
            if row.strand == '+':
                front_overlay.append(front_output[i])
                back_overlay.append(back_output[i])
            else:
                front_overlay.append(back_output[i])
                back_overlay.append(front_output[i])
        overlay_output = np.hstack((front_overlay,overlay,back_overlay))
    else:
        overlay_output = overlay

    return overlay_output