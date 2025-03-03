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
def coordinate_table_to_bed(coordinate_table: pd.DataFrame|str,
                    save_to_file:bool = False,
                    export_as_bedtool:bool = False):
    if isinstance(coordinate_table, str):
        if (os.path.isfile(coordinate_table) == True):
            coordinate_local = pd.read_csv(coordinate_table, sep='\t')
        else:
            logger.error('coordinate_table file not found')
    else:
        coordinate_local = coordinate_table

    coordinate_bed = coordinate_local.iloc[:,0:3]
    columns = coordinate_bed.columns
    new_names = ['chrom', 'start', 'end']
    coordinate_bed.rename(columns={columns[i]: new_names[i] for i in range(3)}, inplace=True)
    coordinate_bed['name'] = coordinate_local.meta_id
    coordinate_bed['score'] = 10
    coordinate_bed['strand'] = coordinate_local.strand
    if save_to_file == True:
        if isinstance(save_to_file, str):
            output_filepath = save_to_file
        else:
            if isinstance(coordinate_table, str):
                output_dir = '/'.join(str.split(coordinate_local, sep ='/')[:-1])
            else:
                output_dir = os.path.dirname(os.path.abspath(__file__)) 
            output_filepath = f'{output_dir}/coordinate.bed'
        coordinate_bed.to_csv(output_filepath, sep='\t', index= False)
    if export_as_bedtool == True:
        coordinate_bed=pybedtools.BedTool.from_dataframe(coordinate_bed)
    return coordinate_bed
#%%
def intersect_process(intersect_df: pd.DataFrame,
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
def bed_io(coordinate_table: pd.DataFrame|str, 
              bed: pd.DataFrame|pybedtools.BedTool|str, 
              save_to_file:bool = False, 
              generate_new_id:bool = False,
              strand_overlap = True,):

    if isinstance(coordinate_table, str):
        if (os.path.isfile(coordinate_table) == True):
            coordinate_local = pd.read_csv(coordinate_table, sep='\t', header=None)
        else:
            logger.error('coordinate_table file not found')
    else:
        coordinate_local = coordinate_table

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
    if generate_new_id == True:
        meta_id = [f'entry_{index}' for index in coordinate_local.index.astype(str)]
        coordinate_local['meta_id'] = meta_id
    else:
        coordinate_local['meta_id'] = coordinate_local.iloc[:,3]
        meta_id = coordinate_local.meta_id.unique()
    

    coordinate_bed=pybedtools.BedTool.from_dataframe(coordinate_local.iloc[:,0:6])
    bed = pybedtools.BedTool(bed)

    intersect_bed=coordinate_bed.intersect(bed_file,loj=True, wa = True, wb =True, s=strand_overlap)
    intersect_df=intersect_bed.to_dataframe()
    with ProcessPoolExecutor(max_workers=40) as executor:
        results  = executor.map(intersect_process, repeat(intersect_df), meta_id)
    bed_out = []
    for result in results:
        bed_out.append(result)
    if save_to_file == True:
        if isinstance(save_to_file, str):
            output_filepath = save_to_file
        else:
            if isinstance(coordinate_table, str):
                output_dir = '/'.join(str.split(coordinate_table, sep ='/')[:-1])
            else:
                output_dir = os.path.dirname(os.path.abspath(__file__))  
            output_filepath = f'{output_dir}/bed_out.p'
        import compress_pickle
        compress_pickle.dump(bed_out, output_filepath, compression="lzma")
        logger.info('done, saving bed_out at: ',output_filepath)
    else:
        logger.info('done, returning bed_out as object')
        return bed_out