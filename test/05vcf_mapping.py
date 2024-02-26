#%%
import sys
import pandas as pd
import numpy as np
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import mapper
from ma_mapper import fetch_data
#%%
input_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/custom_id.fasta.aligned'
aligned_parsed = mapper.parse_alignment(input_filepath, save_to_file= False)
metadata_aligned = mapper.extract_metadata_from_alignment(input_filepath)
#%%
metadata_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11a_coord_with_id.txt'
vcf_mapped=fetch_data.fetch_vcf(metadata_input= metadata_filepath, vcf_input='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/vcf-gnomad/', vcf_format = 'gnomad', custom_id= True)
# %%
metadata_df = pd.read_csv(metadata_filepath, sep='\t')
original_order = metadata_df.iloc[:,4].unique()
#%%
vcf_mapped_sorted = []
for idx, row in metadata_aligned.iterrows():
    vcf_mapped_sorted.append(vcf_mapped[np.where(original_order == row.id)[0][0]])
#%%
filters=mapper.create_filter(aligned_parsed)
# %%
aligned_vcf_overlay=mapper.map_data(vcf_mapped_sorted, aligned_parsed, filters = filters)
#%%
two_spread_color_scale = [[0, "white"],
                    [0.0000000000000000001, "#77ccff"],
                    [0.02, "#55aaff"],
                    [0.04, "#3388ff"],
                    [0.06, "#0066ff"],
                    [0.08, "#0044ff"],
                    [0.1, "#0000ff"],
                    [0.1,"yellow"],
                    [1.0, "red"]]
#%%
import plotly.graph_objects as go
subfamily = 'MER11A'
fig = go.Figure()
fig.add_trace(go.Heatmap(z= aligned_vcf_overlay,
                    colorscale=two_spread_color_scale,
                    colorbar=dict(
                        tick0=0,
                        dtick=1
                       )
))
fig.update_coloraxes(showscale=False)
fig.update_layout(xaxis_showgrid=False, xaxis_zeroline=False,violingap=0,violinmode='overlay',
    title=subfamily + ' vcf overlay',
    xaxis_title="Position (bp)",
    yaxis_title="Sequences",
    legend_title="base",
    showlegend=False,
    )
fig.update_layout(showlegend=False,)
fig.show()
#%%
np.median(aligned_vcf_overlay[aligned_vcf_overlay>0])
# %%
np.std(aligned_vcf_overlay[aligned_vcf_overlay>0])
# %%
