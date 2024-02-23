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
# %%
metadata_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11a_coord_with_id.txt'
maf_mapped=fetch_data.fetch_maf(metadata_input= metadata_filepath, maf_input='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/241genomes/241-mammalian-2020v2b.maf', separated_maf = True,target_species = 'Homo_sapiens', custom_id= True)
#%%
metadata_df = pd.read_csv(metadata_filepath, sep='\t')
original_order = metadata_df.iloc[:,4].unique()
#%%
maf_mapped_sorted = []
for idx, row in metadata_aligned.iterrows():
    maf_mapped_sorted.append(maf_mapped[np.where(original_order == row.id)[0][0]])
#%%
filters=mapper.create_filter(aligned_parsed)
# %%
aligned_maf_overlay=mapper.map_data(maf_mapped_sorted, aligned_parsed, filters = filters)
#%%
import plotly.graph_objects as go
subfamily = 'MER11A'
fig = go.Figure()
fig.add_trace(go.Heatmap(z= aligned_maf_overlay,
                    #colorscale=five_color_scale,
                    #colorbar=dict(
                    #    tick0=0,
                    #    dtick=1
                     #   )
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

# %%
print(np.median(aligned_maf_overlay))
print(np.std(aligned_maf_overlay))
# %%
