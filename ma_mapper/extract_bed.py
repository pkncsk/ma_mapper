#%%
import pybedtools
import numpy as np
import pandas as pd
import os
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat 
#%%
def extract_bed(bed_file, chrom, start_list, end_list, strand):
    drawings = []

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
            bed_file = pybedtools.BedTool(bed_input)
        else:
            print('bed file not found')
    else:
        bed_file = bed_input
    print('extract from bam file: '+ bed_input)
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
        print('done, saving bed_out at: '+output_filepath)
    else:
        print('done, returning bed_out as object')
        return bed_out