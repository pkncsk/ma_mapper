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
bed_mapped=fetch_data.fetch_bed(metadata_input=metadata_filepath, bed_input='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/tfbs_homer/Gata6(Zf).bed', custom_id=True)
#%%
metadata_df = pd.read_csv(metadata_filepath, sep='\t')
original_order = metadata_df.iloc[:,4].unique()
#%%
bed_mapped_sorted = []
for idx, row in metadata_aligned.iterrows():
    bed_mapped_sorted.append(bed_mapped[np.where(original_order == row.id)[0][0]])
#%%
filters=mapper.create_filter(aligned_parsed)
#%%
aligned_bed_overlay=mapper.map_data(bed_mapped_sorted, aligned_parsed, filters = filters)
#%%
import plotly.graph_objects as go
subfamily = 'MER11A'
fig = go.Figure()
fig.add_trace(go.Heatmap(z= aligned_bed_overlay,
                    colorscale="Blues", zauto= False, zmax = 1
                    #colorbar=dict(tick0=0,dtick=1)
))
fig.update_coloraxes(showscale=False)
fig.update_layout(xaxis_showgrid=False, xaxis_zeroline=False,violingap=0,violinmode='overlay',
    title=subfamily + ' AP-1 bed overlay',
    xaxis_title="Position (bp)",
    yaxis_title="Sequences",
    legend_title="base",
    showlegend=False,
    )
fig.update_layout(showlegend=False,)
fig.show()
# %%
bed_filepath ='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/tfbs_homer/Gata6(Zf).bed'
metadata_bed=pybedtools.BedTool.from_dataframe(metadata_df)
bed_file = pybedtools.BedTool(bed_filepath)
# %%
temp_intersect=metadata_bed.intersect(bed_file,wa = True, wb =True, s=True)
# %%
for i in original_order:
    metadata_by_id = metadata_df[metadata_df.id == i]
#%%
metadata_by_id = metadata_df[metadata_df.id == 'n72']
chrom = metadata_by_id.iloc[:,0].unique()[0]
start_list = metadata_by_id.iloc[:,1].to_list()
end_list = metadata_by_id.iloc[:,2].to_list()
strand = metadata_by_id.iloc[:,3].unique()[0]

drawings = []
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
        start_intersect = int(intersect[7])
        end_intersect = int(intersect[8])
        score_intersect = float(intersect[10])
        start_intersect_position = start_intersect - start
        end_intersect_postion = end_intersect - start +1
        canvas[start_intersect_position:end_intersect_postion] += score_intersect
    if strand == '-':
        canvas = np.flip(canvas)
    drawings.append(canvas)
if strand == '+':
    bed_out = np.concatenate(drawings)
else:
    reverse_order=np.flip(drawings, 0)
    bed_out = np.concatenate(reverse_order)
# %%
metadata_by_id_bed=pybedtools.BedTool.from_dataframe(metadata_by_id)
#%%
metadata_bed.intersect(bed_file, wb =True, s=True).head()
#%%
meta_id=metadata_df.id.unique()
metadata_for_bed=metadata_df.iloc[:,0:3]
metadata_for_bed['name'] = metadata_df.id
metadata_for_bed['score'] = 10
metadata_for_bed['strand'] = metadata_df.strand
# %%
bed_filepath ='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/tfbs_homer/Gata4(Zf).bed'
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
                intersect_start_position = intersect_start - canvas_start
                intersect_end_postion = intersect_end - canvas_start +1
                canvas[intersect_start_position:intersect_end_postion] += intersect_score
        if strand == '-':
            canvas = np.flip(canvas)
        drawings.append(canvas)
    if strand == '+':
        bed_out.append(np.concatenate(drawings))
    else:
        reverse_order=np.flip(drawings, 0)
        bed_out.append(np.concatenate(reverse_order))
