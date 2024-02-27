#%%
import sys
import pandas as pd
import numpy as np
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import mapper
from ma_mapper import fetch_sequence
#%%
input_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/custom_id.fasta.aligned'
aligned_parsed = mapper.parse_alignment(input_filepath, save_to_file= False)
metadata_aligned = mapper.extract_metadata_from_alignment(input_filepath)
# %%
metadata_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11a_coord_with_id.txt'
#%%
metadata_df = pd.read_csv(metadata_filepath, sep='\t')
original_order = metadata_df.iloc[:,4].unique()
#%% find border for each id
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
#%% make new metadata 
temp_dict = {'chrom':chrom_list,'start':low_border_list,'end':low_border_list,'strand':strand_list,'id':original_order}
low_border_metadata = pd.DataFrame(temp_dict)
low_border_metadata.start = low_border_metadata.start-500
# %%
source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/hg38.fa'
low_border_records=fetch_sequence.fetch_sequence(low_border_metadata,source_fasta, custom_id= False)

# %%
temp_dict = {'chrom':chrom_list,'start':high_border_list,'end':high_border_list,'strand':strand_list,'id':original_order}
high_border_metadata = pd.DataFrame(temp_dict)
high_border_metadata.end = high_border_metadata.end+500
#%%
source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/hg38.fa'
high_border_records=fetch_sequence.fetch_sequence(high_border_metadata,source_fasta, custom_id= False)
# %%
front_list = []
back_list = []
for idx, strand in enumerate(strand_list):
    if strand == '+':
        front_list.append(low_border_records[idx])
        back_list.append(high_border_records[idx])
    else:
        front_list.append(high_border_records[idx])
        back_list.append(low_border_records[idx])
    
# %%
front_parsed = mapper.parse_alignment(front_list, save_to_file= False)
back_parsed = mapper.parse_alignment(back_list, save_to_file= False)
#%%
filters=mapper.create_filter(aligned_parsed)
row_filter = filters[0]
col_filter = filters[1]
aligned_col_filtered=aligned_parsed[np.ix_(range(len(aligned_parsed)),col_filter)]

# %%
fused_parsed = []
for i in range(len(aligned_parsed)):
    fused_parsed.append(np.concatenate((front_parsed[i], aligned_col_filtered[i], back_parsed[i])))
fused_parsed = np.array(fused_parsed)
#%%
fused_parsed_row_filtered = fused_parsed[np.ix_(row_filter,range(fused_parsed.shape[1]))]
# %%
five_color_scale = [[0, "grey"],
                    [0.2, "grey"],
                    [0.2, "green"],
                    [0.4, "green"],
                    [0.4, "yellow"],
                    [0.6, "yellow"],
                    [0.6, "red"],
                    [0.8, "red"],
                    [0.8, "blue"],
                    [1.0, "blue"]]
six_color_scale = [[0, "grey"],
                            [0.16666666666666666, "grey"],

                            [0.16666666666666666, "green"],
                            [0.16666666666666666*2, "green"],

                            [0.16666666666666666*2, "yellow"],
                            [0.16666666666666666*3, "yellow"],

                            [0.16666666666666666*3, "red"],
                            [0.16666666666666666*4, "red"],

                            [0.16666666666666666*4, "blue"],
                            [0.16666666666666666*5, "blue"],

                            [0.16666666666666666*5, "black"],
                            [1.0, "black"]]
#%%
import plotly.graph_objects as go
subfamily = 'MER11A'
fig = go.Figure()
fig.add_trace(go.Heatmap(z= fused_parsed_row_filtered,
                    colorscale=five_color_scale,
                    colorbar=dict(
                        tick0=0,
                        dtick=1
                        )
))
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
#%%