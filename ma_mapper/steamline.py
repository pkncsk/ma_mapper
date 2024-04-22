import sys
import pandas as pd
import numpy as np
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import mapper
from ma_mapper import fetch_data
from ma_mapper import fetch_sequence
#%%
def metadata_extend(metadata_df, extension_length = 500, export_strand_list = False):
    original_order = metadata_df.iloc[:,4].unique()
    chrom_list = []
    low_border_list = []
    high_border_list = []
    strand_list = []
    for uniq_meta_id in original_order:
        metadata_by_id = metadata_df[metadata_df.id == uniq_meta_id]
        chrom_list.append(metadata_by_id.iloc[:,0].unique()[0])
        low_border_list.append(min(metadata_by_id.iloc[:,1]))
        high_border_list.append(max(metadata_by_id.iloc[:,2]))
        strand_list.append(metadata_by_id.iloc[:,3].unique()[0])
    
    temp_dict = {'chrom':chrom_list,'start':low_border_list,'end':low_border_list,'strand':strand_list,'id':original_order}
    low_border_metadata = pd.DataFrame(temp_dict)
    low_border_metadata.start = low_border_metadata.start-extension_length

    temp_dict = {'chrom':chrom_list,'start':high_border_list,'end':high_border_list,'strand':strand_list,'id':original_order}
    high_border_metadata = pd.DataFrame(temp_dict)
    high_border_metadata.end = high_border_metadata.end+extension_length
    if export_strand_list == False:
        return low_border_metadata, high_border_metadata
    else:
        return low_border_metadata, high_border_metadata, strand_list

def alignment_matrix(alignment_filepath, metadata_filepath, source_fasta, extension_length = 500):
    aligned_parsed = mapper.parse_alignment(alignment_filepath, save_to_file= False)
    metadata_aligned = mapper.extract_metadata_from_alignment(alignment_filepath)
    metadata_df = pd.read_csv(metadata_filepath, sep='\t')
    original_order = metadata_df.iloc[:,4].unique()
    low_border_metadata, high_border_metadata, strand_list= metadata_extend(metadata_df, extension_length, export_strand_list=True)

    low_border_records=fetch_sequence.fetch_sequence(low_border_metadata,source_fasta, custom_id= False)
    high_border_records=fetch_sequence.fetch_sequence(high_border_metadata,source_fasta, custom_id= False)

    front_list = []
    back_list = []
    for idx, strand in enumerate(strand_list):
        if strand == '+':
            front_list.append(low_border_records[idx])
            back_list.append(high_border_records[idx])
        else:
            front_list.append(high_border_records[idx])
            back_list.append(low_border_records[idx])
    front_parsed = mapper.parse_alignment(front_list, save_to_file= False)
    back_parsed = mapper.parse_alignment(back_list, save_to_file= False)
    filters=mapper.create_filter(aligned_parsed,row_threshold=0.1, col_threshold=0.1,col_content_threshold=0.1)
    filters=mapper.create_filter(aligned_parsed,row_threshold=0.1, col_threshold=0.1,col_content_threshold=0.1)
    row_filter = filters[0]
    col_filter = filters[1]
    aligned_col_filtered=aligned_parsed[np.ix_(range(len(aligned_parsed)),col_filter)]

    front_parsed_sorted = []
    back_parsed_sorted = []
    for idx, row in metadata_aligned.iterrows():
        front_parsed_sorted.append(front_parsed[np.where(original_order == row.id)[0][0]])
        back_parsed_sorted.append(back_parsed[np.where(original_order == row.id)[0][0]])
    fused_parsed = []
    for i in range(len(aligned_parsed)):
        fused_parsed.append(np.concatenate((front_parsed_sorted[i], aligned_col_filtered[i], back_parsed_sorted[i])))
    fused_parsed = np.array(fused_parsed)

    fused_parsed_row_filtered = fused_parsed[np.ix_(row_filter,range(fused_parsed.shape[1]))]
    metadata_aligned_filtered=metadata_aligned.iloc[row_filter,:] 
    return fused_parsed_row_filtered, metadata_aligned_filtered    

def bigwig_matrix(alignment_filepath, metadata_filepath,bigwig_filepath, extension = False, extension_length = 500):
    aligned_parsed = mapper.parse_alignment(alignment_filepath, save_to_file= False)
    metadata_aligned = mapper.extract_metadata_from_alignment(alignment_filepath)
    metadata_df = pd.read_csv(metadata_filepath, sep='\t')
    original_order = metadata_df.iloc[:,4].unique()
    if extension:
        low_border_metadata, high_border_metadata, strand_list= metadata_extend(metadata_df, extension_length, export_strand_list=True)
    phylop_mapped=fetch_data.fetch_bigwig(metadata_input= metadata_filepath, bigwig_input=bigwig_filepath, custom_id=True)

    bigwig_mapped_sorted = []
    for idx, row in metadata_aligned.iterrows():
        bigwig_mapped_sorted.append(phylop_mapped[np.where(original_order == row.id)[0][0]])
    
    filters=mapper.create_filter(aligned_parsed)
    row_filter = filters[0]
    col_filter = filters[1]
    aligned_filtered=aligned_parsed[np.ix_(row_filter,col_filter)]
    aligned_bigwig_overlay=mapper.map_data(bigwig_mapped_sorted, aligned_parsed, filters = filters)
    metadata_aligned_filtered=metadata_aligned.iloc[row_filter,:]
    
    if extension == False:
        return aligned_bigwig_overlay, metadata_aligned_filtered
    else:
        low_border_bigwig_mapped=fetch_data.fetch_bigwig(metadata_input= low_border_metadata, bigwig_input=bigwig_filepath, custom_id= True)
        high_border_bigwig_mapped=fetch_data.fetch_bigwig(metadata_input= high_border_metadata, bigwig_input=bigwig_filepath, custom_id= True)
        bigwig_front_list = []
        bigwig_back_list = []
        for idx, strand in enumerate(strand_list):
            if strand == '+':
                bigwig_front_list.append(low_border_bigwig_mapped[idx])
                bigwig_back_list.append(high_border_bigwig_mapped[idx])
            else:
                bigwig_front_list.append(high_border_bigwig_mapped[idx])
                bigwig_back_list.append(low_border_bigwig_mapped[idx])
        bigwig_front_sorted = []
        bigwig_back_sorted = []
        for idx, row in metadata_aligned.iterrows():
            if row_filter[idx]:
                bigwig_front_sorted.append(bigwig_front_list[np.where(original_order == row.id)[0][0]])
                bigwig_back_sorted.append(bigwig_back_list[np.where(original_order == row.id)[0][0]])
        fused_bigwig_mapped = []
        for i in range(len(aligned_bigwig_overlay)):
            fused_bigwig_mapped.append(np.concatenate((bigwig_front_sorted[i], aligned_bigwig_overlay[i], bigwig_back_sorted[i])))
        fused_bigwig_mapped = np.array(fused_bigwig_mapped)

        return fused_bigwig_mapped, metadata_aligned_filtered

def matrix_cluster(matrix, extension = False, extension_length = 500, export_order = False, clustering_metric = 'euclidean', clustering_method ='ward'):
    import scipy
    if extension == True:
        guide_matrix = matrix[:,extension_length:-extension_length]
    martix_treated=np.nan_to_num(guide_matrix)
    distance_matrix = scipy.spatial.distance.pdist(martix_treated, metric=clustering_metric)
    linkage_array = scipy.cluster.hierarchy.linkage(distance_matrix, method = clustering_method)
    optimised_order=scipy.cluster.hierarchy.optimal_leaf_ordering(linkage_array, distance_matrix, metric=clustering_metric)
    new_order = scipy.cluster.hierarchy.leaves_list(optimised_order)
    matrix_sorted = []
    for i in new_order:
        matrix_sorted.append(matrix[i])
    clustered_matrix = np.array(matrix_sorted)
    if export_order == True:
        return clustered_matrix, new_order
    else:
        return clustered_matrix

def normalise(alignment_df, target_df):
    normalisation_mask = np.count_nonzero(alignment_df, axis=0)
    target_df_treated=np.nan_to_num(target_df)
    sum_target_df_treated = np.sum(target_df_treated, axis = 0)
    target_df_averaged = sum_target_df_treated/normalisation_mask
    return target_df_averaged