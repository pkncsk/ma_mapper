#%%
import pybedtools
import numpy as np
import pandas as pd
import os
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat 
from . import logger
import sys
if sys.version_info >= (3, 8, 0):
    from typing import Literal, Tuple, List
else:
    from typing_extensions import Literal, Tuple, List
#%%
def metadata_to_bed(metadata: pd.DataFrame|str,
                    output_dir: str|None = None, 
                    save_to_file:bool = False,
                    export_bedtool:bool = False):
    if isinstance(metadata, str):
        if (os.path.isfile(metadata) == True):
            metadata_local = pd.read_csv(metadata, sep='\t')
        else:
            logger.error('metadata file not found')
    else:
        metadata_local = metadata
    if output_dir is None:
        if isinstance(metadata, str):
            output_dir = '/'.join(str.split(metadata_local, sep ='/')[:-1])
        else:
            output_dir = os.path.dirname(os.path.abspath(__file__))  
    metadata_bed = metadata_local.iloc[:,0:3]
    columns = metadata_bed.columns
    new_names = ['chrom', 'start', 'end']
    metadata_bed.rename(columns={columns[i]: new_names[i] for i in range(3)}, inplace=True)
    metadata_bed['name'] = metadata_local.meta_id
    metadata_bed['score'] = 10
    metadata_bed['strand'] = metadata_local.strand
    if save_to_file == True:
        metadata_bed.to_csv(output_dir, sep='\t', index= False)
    if export_bedtool:
        metadata_bed=pybedtools.BedTool.from_dataframe(metadata_bed)
    return metadata_bed

#%%

def _extract_bed(bed_file, chrom, start_list, end_list, strand):
    drawings = []

    if isinstance(bed_file, str):
        if (os.path.isfile(bed_file) == True):
            bed_file_local = pybedtools.BedTool(bed_file)
        else:
            logger.error('bed file not found')
    else:
        bed_file_local = bed_file
    for i in range(len(start_list)):
        start = start_list[i]
        end = end_list[i]
        prep_dict = {'chrom':[chrom], 'start':[start], 'end':[end],'name':['temp'],'score':[10], strand:[strand]}
        temp_df = pd.DataFrame(prep_dict)
        metadata_bed=pybedtools.BedTool.from_dataframe(temp_df)
        temp_intersect=metadata_bed.intersect(bed_file_local,wa=True, wb =True, s=True)
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
def extract_bed(intersect_df: pd.DataFrame,
                unique_id):
    intersect_df_local = intersect_df
    intersect_by_id=intersect_df_local[intersect_df_local['name'] == unique_id]
    unique_start=intersect_by_id.start.unique()
    strand=intersect_by_id.strand.unique()[0]
    drawings = []
    for start in unique_start:
        intersect_by_id_by_start=intersect_by_id[intersect_by_id.start == start]
        canvas_start = intersect_by_id_by_start.start.unique()[0]
        canvas_end = intersect_by_id_by_start.end.unique()[0]
        if(len(intersect_by_id_by_start.end.unique())>1):
            logger.info('multiple length frags')
        canvas = np.zeros(canvas_end-canvas_start)
        for idx,row in intersect_by_id_by_start.iterrows():
            if int(row.thickEnd) != -1:
                intersect_start = int(row.thickEnd)
                intersect_end = int(row.itemRgb)
                intersect_score = int(row.blockSizes)
                #capping out of bound position
                if intersect_end > canvas_end:
                    intersect_end = canvas_end
                if intersect_start < canvas_start:
                    intersect_start = canvas_start
                if strand == '+':
                    start_intersect_position = intersect_start - canvas_start
                    end_intersect_postion = intersect_end - canvas_start

                else:
                    start_intersect_position = canvas_end - intersect_end
                    end_intersect_postion = canvas_end - intersect_start
                canvas[start_intersect_position:end_intersect_postion] += intersect_score
        #if strand == '-':
            #canvas = np.flip(canvas)
        drawings.append(canvas)
    if strand == '-':
        drawings.reverse()
    bed_out=np.concatenate(drawings)
    return bed_out
#%%
_MODE = Literal['legacy','dirty']
def bed_io(metadata: pd.DataFrame|str, 
              bed: pd.DataFrame|pybedtools.BedTool|str, 
              output_dir: str|None = None, 
              save_to_file:bool = False, 
              custom_id:bool = False,
              custom_prefix = 'entry',
              strand_overlap = True,
              mode:_MODE = 'dirty'):

    if isinstance(metadata, str):
        if (os.path.isfile(metadata) == True):
            metadata_local = pd.read_csv(metadata, sep='\t', header=None)
        else:
            logger.error('metadata file not found')
    else:
        metadata_local = metadata
    if output_dir is None:
        if isinstance(metadata, str):
            output_dir = '/'.join(str.split(metadata, sep ='/')[:-1])
        else:
            output_dir = os.path.dirname(os.path.abspath(__file__))  
    if isinstance(bed, str):
        if (os.path.isfile(bed) == True):
            bed_file = pybedtools.BedTool(bed)
            logger.info('extract from bed file: '+ bed)
        else:
            logger.error('bed file not found')
    elif isinstance(bed, pd.DataFrame):
        bed_file=pybedtools.BedTool.from_dataframe(bed.iloc[:,0:6])
        logger.info('extract from dataframe')
    else:
        bed_file = bed
    if custom_id == False:
        meta_id = [f'{custom_prefix}_{index}' for index in metadata_local.index.astype(str)]
        metadata_local['meta_id'] = meta_id
    else:
        metadata_local['meta_id'] = metadata_local.iloc[:,3]
        meta_id = metadata_local.meta_id.unique()
    
    
    if mode == 'legacy':
        grouped = metadata_local.groupby('meta_id', sort=False)
        chrom_list = grouped.apply(lambda x: x.iloc[:,0].unique()[0], include_groups=False).tolist()
        start_list = grouped.apply(lambda x: x.iloc[:,1].tolist(), include_groups=False).tolist()
        end_list = grouped.apply(lambda x: x.iloc[:,2].tolist(), include_groups=False).tolist()
        strand_list = grouped.apply(lambda x: x.iloc[:,5].unique()[0], include_groups=False).tolist()
        with ProcessPoolExecutor(max_workers=40) as executor:
            results  = executor.map(_extract_bed, repeat(bed_file), chrom_list, start_list,end_list,strand_list)
    elif mode == 'dirty':
        metadata_bed=pybedtools.BedTool.from_dataframe(metadata_local.iloc[:,0:6])
        bed = pybedtools.BedTool(bed)

        intersect_bed=metadata_bed.intersect(bed_file,loj=True, wa = True, wb =True, s=strand_overlap)
        intersect_df=intersect_bed.to_dataframe()
        with ProcessPoolExecutor(max_workers=40) as executor:
            results  = executor.map(extract_bed, repeat(intersect_df), meta_id)
    bed_out = []
    for result in results:
        bed_out.append(result)
    if save_to_file == True:
        import compress_pickle
        output_filepath = f'{output_dir}/bed_out.lzma'
        compress_pickle.dump(bed_out, output_filepath, compression="lzma")
        logger.info('done, saving bed_out at: ',output_filepath)
    else:
        logger.info('done, returning bed_out as object')
        return bed_out