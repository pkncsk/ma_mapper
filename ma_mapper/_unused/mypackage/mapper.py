#%%
import pandas as pd
import numpy as np
import sys
from . import logger
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
            output_filepath = os.path.dirname(os.path.abspath(__file__)) + '/metadata_aligned.txt'
        new_metadata.to_csv(output_filepath, sep='\t', index= False)
        print('save new metadata at', output_filepath)
    else:
        return new_metadata

def parse_alignment(alignment_file, save_to_file=False, output_file=None):
    if output_file is None:
        if isinstance(alignment_file, str):
            output_file = alignment_file + '.parsed'
        else:
            output_file = os.path.dirname(os.path.abspath(__file__)) +'/'+ 'alignment.parsed'
    output_dir = '/'.join(str.split(output_file, sep ='/')[:-1])
    logger.info('parse alignment')
    if isinstance(alignment_file, str):
        if (os.path.isfile(alignment_file) == True):
            records = load_alignment_file(alignment_file)
        else:
            print('alignment file not found')
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
            print('alignment file not found')
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
        logger.info('done, saving parsed array at: '+output_file)
    else: 
        return result_array
#%%
def import_parsed_array(parsed_array_file):
    import h5py
    with h5py.File(import_parsed_array+'/parsed_array.h5', 'r') as hf:
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
def map_data(data_file, sorted_parsed_array, filters= None):
    if isinstance(data_file, str) == True:
        import compress_pickle
        data_file = compress_pickle.load(data_file)
    if filters is not None:
        if isinstance(filters,str) == True:
            import pickle
            filters = pickle.load(open(filters)) 
        else:
            filters = filters

    canvas = np.zeros(sorted_parsed_array.shape, np.float32)
    for idx, row in enumerate(sorted_parsed_array):
        canvas[idx, row>0] = data_file[idx]
    mapped_data = canvas
    if filters is not None:
        row_filter = filters[0]
        col_filter = filters[1]
        mapped_data=mapped_data[np.ix_(row_filter,col_filter)]
    return mapped_data
#%%
#FIXME: wrong nonzero
#from typing import Literal <- python3.8+
if sys.version_info >= (3, 8, 0):
    from typing import Literal
else:
    from typing_extensions import Literal
_METHOD = Literal['average',]
_MODE = Literal['all','present']
def normalise(mapped_data, method:_METHOD = 'average', mode:_MODE = 'present', outlier = None):
    #if isinstance(mapped_data, str) == True:
    if outlier is not None:
        mapped_data[mapped_data >= outlier] = 0
    if method == 'average':
        numerator =np.sum(mapped_data, axis = 0)
        if mode == 'present':
            denominator = np.count_nonzero(mapped_data, axis=0)
        else:
            denominator = np.shape(mapped_data)[0] 
    normalised_array = numerator/denominator
    return normalised_array
#steamline functions
def parse_and_filter(alignment_file):
    aligned_parsed = parse_alignment(alignment_file, save_to_file= False)
    metadata_aligned = extract_metadata_from_alignment(alignment_file)
    metadata_aligned['original_order'] = metadata_aligned.index
    filters=create_filter(aligned_parsed)
    row_filter = filters[0]
    col_filter = filters[1]
    aligned_filtered=aligned_parsed[np.ix_(row_filter,col_filter)]
    metadata_aligned_filtered=metadata_aligned.iloc[row_filter,:]
    return aligned_filtered, metadata_aligned_filtered

def match_age_to_id_metadata(metadata_df):
    species_reference = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/combined_age_div/combined_age_and_div.txt'
    main_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
    age_div_table = pd.read_csv(species_reference, sep='\t')
    age_div_table['id'] = age_div_table.repName+'_'+age_div_table.internal_id.astype(str)
    te_age=age_div_table[['id','te_age','te_div']].drop_duplicates()
    metadata_with_te_age=metadata_df.merge(te_age, on = 'id', how ='left')
    return metadata_with_te_age