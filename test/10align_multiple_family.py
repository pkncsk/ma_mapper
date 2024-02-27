#%%
from curses import meta
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
#%% from age_div table
subfamily = ['MER11A','MER11B','MER11C']
input_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/combined_age_div/combined_age_and_div.txt'
main_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
age_div_table = pd.read_csv(input_filepath, sep='\t')
subfam_table=age_div_table[age_div_table.repName.isin(subfamily)]
subfam_table = subfam_table[subfam_table.genoName.isin(main_chr)]
subfam_coord = subfam_table[['genoName','genoStart','genoEnd','strand']]
#subfam_coord['id'] = 'n'+subfam_table.internal_id.astype(str)
subfam_coord['id'] = subfam_table.repName+'_'+subfam_table.internal_id.astype(str)
#subfam_coord['id'] = 'n'+subfam_table.internal_id.astype(str) + '_' + subfam_table.index.astype(str)
# %%
subfam_coord.to_csv('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11_coord_with_id.txt', sep='\t', index= False)
#%%

import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import fetch_sequence
from ma_mapper import mafft_align
from ma_mapper import fetch_data
#%%

source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_fasta/hg38.fa'
metadata = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11_coord_with_id.txt'
fasta_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11.fasta'

fetch_sequence.fetch_sequence(metadata,source_fasta,output_filepath =fasta_file, save_to_file= True,custom_id= True)
#%%
mafft_align.mafft_wrapper(fasta_file)
# %%
from ma_mapper import mapper
import numpy as np
import pandas as pd
#%%
input_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11.fasta.aligned'
aligned_parsed = mapper.parse_alignment(input_filepath, save_to_file= False)
metadata_aligned = mapper.extract_metadata_from_alignment(input_filepath)
#%%
filters=mapper.create_filter(aligned_parsed)
# %%
row_filter = filters[0]
col_filter = filters[1]
aligned_filtered=aligned_parsed[np.ix_(row_filter,col_filter)]
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
fig.add_trace(go.Heatmap(z= aligned_filtered,
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