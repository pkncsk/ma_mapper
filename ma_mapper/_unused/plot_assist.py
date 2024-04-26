#%%
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
from importlib import metadata
import sys
import numpy as np
sys.path.append('/home/pc575/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import mapper
#%%
input_filepath = 'data/_ma_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/custom_id.fasta.aligned'
def plot_alignment(input_filepath, custom_title = None, save_to_file = False, output_file = None):
    aligned_parsed = mapper.parse_alignment(input_filepath, save_to_file= False)
    filters=mapper.create_filter(aligned_parsed)
    row_filter = filters[0]
    col_filter = filters[1]
    aligned_filtered=aligned_parsed[np.ix_(row_filter,col_filter)]
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
    fig = go.Figure()
    fig.add_trace(go.Heatmap(z= aligned_filtered,
                        colorscale=five_color_scale,
                        colorbar=dict(
                            tick0=0,
                            dtick=1
                            )
    ))
    if custom_title is None:
        plot_title = 'Alignment'
    else:
        plot_title = custom_title
    fig.update_coloraxes(showscale=False)
    fig.update_layout(xaxis_showgrid=False, xaxis_zeroline=False,violingap=0,violinmode='overlay',
        title=plot_title,
        xaxis_title="Position (bp)",
        yaxis_title="Sequences",
        legend_title="base",
        showlegend=False,
        )
    fig.update_layout(showlegend=False,)
    if save_to_file == True:
        if output_file is None:
            output_filepath = './align_img.png'
        else:
            output_filepath = output_file
        fig.write_image(output_filepath)
    else:
        fig.show()
        