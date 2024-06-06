#%%
from turtle import back
import pandas as pd
from . import sequence_alignment
import numpy as np
import sys
import os
from . import logger
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
    dump_list = []
    chrom_list = []
    start_list = []
    end_list = []
    strand_list = []
    for record in records:
        metadata_list=record.name.split(sep='::')
        for i in range(len(metadata_list)):
            if metadata_list[i].startswith('chr') == False:
                dump_list.append(metadata_list[i])
            else:
                chrom_list.append(metadata_list[i])
                start_list.append(int(metadata_list[i+1]))
                end_list.append(int(metadata_list[i+2]))
                strand_list.append(metadata_list[i+3])
                break
    new_metadata_dict = {'chrom':chrom_list,'start':start_list,'end':end_list,'strand':strand_list,'id':dump_list}
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
def map_data(data_file:str,
             sorted_parsed_array: np.ndarray, 
             filters: bool = True,
             custom_filters = None,
             col_threshold: int = 0.50, 
             col_content_threshold:int = 0.10, 
             row_threshold:int = 0.50,)-> np.ndarray:
    if isinstance(data_file, str) == True:
        import compress_pickle
        data_file = compress_pickle.load(data_file)
    def check_list_type(data, type_check):
        return all(isinstance(item, type_check) for sublist in data for item in sublist)
    if check_list_type(sorted_parsed_array, (int, float, type(np.nan))):
        canvas = np.zeros(sorted_parsed_array.shape, np.float32)
    else:
        canvas = np.zeros(sorted_parsed_array.shape, dtype=object)
    
    for idx, row in enumerate(sorted_parsed_array):
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
_NORM_METHOD = Literal['average',]
def normalise(alignment: str|pd.DataFrame,
              mapped_data: pd.DataFrame, 
              method:_NORM_METHOD = 'average', 
              ):
    if method == 'average':
        normalizer = np.count_nonzero(alignment, axis=0)
        data_sum = np.nansum(mapped_data, axis=0)
        if method == 'average':
            normalised_data = data_sum/normalizer 
    
    return normalised_data

def flank_sequence_io(metadata: pd.DataFrame,
                    source_fasta: str|None = None,
                    metadata_out:bool = False, 
                    extension_length:int = 500) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if source_fasta is None:
        source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/hg38.fa'

    front_metadata = metadata.reset_index(drop=True).copy()
    front_metadata.end = front_metadata.start
    front_metadata.start = front_metadata.start - extension_length
    
    back_metadata = metadata.reset_index(drop=True).copy()
    back_metadata.start = back_metadata.end 
    back_metadata.end = back_metadata.end + extension_length
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
        cf_kwargs = {key: value for key, value in kwargs.items() if key in ('col_threshold','col_content_threshold', 'row_threshold')}
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
                             reference_table:str|None = None) -> pd.DataFrame:
    if isinstance(metadata, str):
        if os.path.isfile(metadata):
            metadata_df = pd.read_csv(metadata, sep='\t')
        else:
            logger.error('metadata file not found')
    elif isinstance(metadata, pd.DataFrame):
        metadata_df = metadata
    if reference_table is None:
        reference_table = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/combined_age_div/combined_age_and_div.txt'
    age_div_table = pd.read_csv(reference_table, sep='\t')
    age_div_table['id'] = age_div_table.repName+'_'+age_div_table.internal_id.astype(str)
    te_age=age_div_table[['id','te_age','te_div']].drop_duplicates()
    metadata_with_te_age=metadata_df.merge(te_age, on = 'id', how ='left')
    return metadata_with_te_age

_FORMAT = Literal['bam_max','bam_min','bam_forward','bam_reverse','bigwig','bed']
def map_and_overlay(aligment:str,
            metadata:str,
            data_file:str,
            data_format:_FORMAT,
            filters:bool = True,
            extension_length: int|None = None,
            *args,**kwargs,
            )->Tuple[np.ndarray,pd.DataFrame,np.ndarray,np.ndarray,np.ndarray,np.ndarray]:

    alignment_parsed, metadata_aligned, filters  = parse_and_filter(alignment_file=aligment, preprocess_out=True, **kwargs)
    
    if data_format in ['read_max','read_min','read_forward','read_reverse', 'normal_max','normal_min','normal_forward','normal_reverse']:
        from . import extract_bam
        output=extract_bam.bam_io(metadata= metadata, bam_file=data_file,bam_format=data_format, custom_id=True,**kwargs)
    elif data_format in ['bigwig']:
        from . import extract_bigwig
        output=extract_bigwig.bigwig_io(metadata=metadata, bigwig=data_file, custom_id=True,**kwargs)
    elif data_format in ['bed']:
        from . import extract_bed
        output=extract_bed.bed_io(metadata=metadata, bed=data_file, custom_id=True,**kwargs)
    elif data_format in ['maf']:
        from . import extract_maf
        output=extract_maf.maf_io(metadata=metadata, maf=data_file,custom_id=True,**kwargs)
    else:
        logger.error('dataformat field unspecified')
    metadata_df = pd.read_csv(metadata, sep='\t')
    source_order = metadata_df.iloc[:,4].unique()
    output_sorted = []
    for idx, row in metadata_aligned.iterrows():
        output_sorted.append(output[np.where(source_order == row.id)[0][0]])
    
    md_kwargs = {key: value for key, value in kwargs.items() if key in ('col_threshold','col_content_threshold', 'row_threshold')} 
    overlay = map_data(data_file=output_sorted, sorted_parsed_array = alignment_parsed, custom_filters=filters, **md_kwargs)
    
    row_filter = filters[0]
    col_filter = filters[1]
    metadata_filtered=metadata_aligned.iloc[row_filter,:]

    if extension_length is not None:
        ffs_kwargs = {key: value for key, value in kwargs.items() if key in ('source_fasta')} 
        front_metadata, back_metadata = flank_sequence_io(metadata_filtered, metadata_out=True, extension_length=extension_length, **ffs_kwargs)
        if data_format in ['bam_max','bam_min','bam_forward','bam_reverse']:
            
            front_output=extract_bam.bam_io(metadata= front_metadata, bam_file=data_file,bam_format=data_format, **kwargs)
            back_output=extract_bam.bam_io(metadata= back_metadata, bam_file=data_file,bam_format=data_format, **kwargs)
        elif data_format in ['bigwig']:
            front_output=extract_bigwig.bigwig_io(metadata=front_metadata, bigwig=data_file,**kwargs)
            back_output=extract_bigwig.bigwig_io(metadata=back_metadata, bigwig=data_file,**kwargs)
        elif data_format in ['bed']:
            front_output=extract_bed.bed_io(metadata=front_metadata, bed=data_file,**kwargs)
            back_output=extract_bed.bed_io(metadata=back_metadata, bed=data_file,**kwargs)
        elif data_format in ['maf']:
            front_output=extract_maf.maf_io(metadata=front_metadata, maf=data_file,custom_id=True,**kwargs)
            back_output=extract_maf.maf_io(metadata=back_metadata, bed=data_file,**kwargs)
        front_overlay = []
        back_overlay = []
        for i, row in metadata_filtered.iterrows():
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