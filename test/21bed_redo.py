#%%
import sys
from tracemalloc import start
import pandas as pd
import numpy as np
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import mapper
from ma_mapper import fetch_data
import pybedtools
#%%
input_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/custom_id.fasta.aligned'
aligned_parsed = mapper.parse_alignment(input_filepath, save_to_file= False)
metadata_aligned = mapper.extract_metadata_from_alignment(input_filepath)
#%%
metadata_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11a_coord_with_id.txt'
#%%
metadata_df = pd.read_csv(metadata_filepath, sep='\t')
original_order = metadata_df.iloc[:,4].unique()
#%%
meta_id=metadata_df.id.unique()
metadata_for_bed=metadata_df.iloc[:,0:3]
metadata_for_bed['name'] = metadata_df.id
metadata_for_bed['score'] = 10
metadata_for_bed['strand'] = metadata_df.strand
# %%
bed_filepath ='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/tfbs_homer/Gata6(Zf).bed'
metadata_bed=pybedtools.BedTool.from_dataframe(metadata_for_bed)
bed_file = pybedtools.BedTool(bed_filepath)
# %%
temp_intersect=metadata_bed.intersect(bed_file,loj=True, wa = True, wb =True, s=True)
intersect_df=temp_intersect.to_dataframe()
#%%
bed_out = []
for unique_id in meta_id:
    intersect_by_id=intersect_df[intersect_df.name == unique_id]
    unique_start=intersect_by_id.start.unique()
    strand=intersect_by_id.strand.unique()[0]
    drawings = []
    for start in unique_start:
        intersect_by_id_by_start=intersect_by_id[intersect_by_id.start == start]
        canvas_start = intersect_by_id_by_start.start.unique()[0]
        canvas_end = intersect_by_id_by_start.end.unique()[0]
        if(len(intersect_by_id_by_start.end.unique())>1):
            print('multiple length frags')
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
    if strand == '+':
        bed_out.append(np.concatenate(drawings))
    else:
        reverse_order=np.flip(drawings, 0)
        bed_out.append(np.concatenate(reverse_order))
#%%
bed_mapped_sorted = []
for idx, row in metadata_aligned.iterrows():
    bed_mapped_sorted.append(bed_out[np.where(original_order == row.id)[0][0]])
#%%
filters=mapper.create_filter(aligned_parsed)
#%%
aligned_bed_overlay=mapper.map_data(bed_mapped_sorted, aligned_parsed, filters = filters)
# %%
row_filter = filters[0]
col_filter = filters[1]
aligned_filtered=aligned_parsed[np.ix_(row_filter,col_filter)]
# %%
text_array=np.array(aligned_filtered).astype('object')
text_array[text_array==0] = ''
text_array[text_array==1] = 'A'
text_array[text_array==2] = 'C'
text_array[text_array==3] = 'T'
text_array[text_array==4] = 'G'
#%%
blue_color_scale = [[0, 'rgba(0, 0, 0, 0.2)'],[1.0, 'rgba(8, 54, 116, 0.5)']]
five_color_scale = [[0, 'rgba(128, 128, 128, 0.2)'],[0.2, 'rgba(128, 128, 128, 0.2)'],[0.2, 'rgba(0, 128, 0, 0.2)'],[0.4, 'rgba(0, 128, 0, 0.2)'],[0.4, 'rgba(255, 255, 0, 0.2)'],[0.6, 'rgba(255, 255, 0, 0.2)'],[0.6, 'rgba(0, 0, 255, 0.2)'],[0.8, 'rgba(0, 0, 255, 0.2)'],[0.8, 'rgba(255, 0, 0, 0.2)'],[1.0, 'rgba(255, 0, 0, 0.2)']]
#%%
#%%
import plotly.graph_objects as go
subfamily = 'MER11A'
fig = go.Figure()


fig.add_trace(go.Heatmap(z= aligned_filtered,
                    colorscale=five_color_scale,
                    text = text_array#colorbar=dict(ick0=0,dtick=1)
))
fig.add_trace(go.Heatmap(z= aligned_bed_overlay,
                    colorscale=blue_color_scale, zauto= False, zmax = 1
                    #colorbar=dict(tick0=0,dtick=1)
))
fig.update_traces(text=text_array, texttemplate="%{text}")
fig.update_coloraxes(showscale=False)
fig.update_layout(xaxis_showgrid=False, xaxis_zeroline=False,violingap=0,violinmode='overlay',
    title=subfamily + ' Alignment',
    xaxis_title="Position (bp)",
    yaxis_title="Sequences",
    legend_title="base",
    showlegend=False,
    )
fig.update_layout(showlegend=False,)
fig.show()

# %%
overlap_id=intersect_df[intersect_df.blockSizes != -1].name.to_list()
# %%
overlap_index=metadata_aligned[metadata_aligned.id.isin(overlap_id)].index
nonoverlap_index=metadata_aligned[~metadata_aligned.id.isin(overlap_id)].index
# %%
