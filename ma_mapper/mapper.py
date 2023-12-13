#%%
import numpy as np
import sys
#%%
def parse_alignment(alignment_file, save_to_file=False, output_file=None):
    if output_file is None:
        output_file = alignment_file + '.parsed'
    import logging
    log_path = output_file+'.parsed.log'
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
    logging.info('parse alignment: '+ alignment_file)
    if isinstance(alignment_file, str) == True:
        from Bio import SeqIO
        records = (r for r in SeqIO.parse(alignment_file, "fasta"))
    else:
        records = alignment_file
    #counting sequence records
    seq_count = 0
    for i, value in enumerate(records):
        if seq_count == 0:
            seq_length = len(value)
        seq_count = seq_count + 1
    if seq_count == 0:
        logging.info('no records')
    else:
        seq_id_list = list()
        for idx, content in enumerate(records):
            seq_id_list.append(content.name)
    #create parsed array
    parsed_array = np.zeros((seq_count, seq_length), dtype=np.uint8)
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
        logging.info('done, saving parsed array at: '+output_file)
    else: 
        return result_array
#%%
def import_parsed_array(parsed_array_file):
    import h5py
    with h5py.File(import_parsed_array+'/parsed_array.h5', 'r') as hf:
        parsed_array = hf['result_array'][:]  
    return parsed_array
#%%
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
def map_data(data_file, sorted_parsed_array, filter= None):
    if isinstance(data_file, str) == True:
        import compress_pickle
        data_file = compress_pickle.load(data_file)
    if filter is not None:
        if isinstance(filter,str) == True:
            import pickle
            filter = pickle.load(open(filter)) 
        else:
            filter = filter

    canvas = np.zeros(sorted_parsed_array.shape, np.float32)
    for index, row in enumerate(sorted_parsed_array):
        canvas[index, row>0] = data_file[index]
    mapped_data = canvas
    if filter is not None:
        row_filter = filter[0]
        col_filter = filter[1]
        mapped_data=mapped_data[np.ix_(row_filter,col_filter)]
    return mapped_data