#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import extract_bigwig
from ma_mapper import mapper
#%%
subfamily='protein_coding_sequences_sliced'
source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_fasta/hg38.fa'
alignment_file = f'/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_main/alignment/{subfamily}_.fasta'
alignment_extended, metadata_filtered= mapper.parse_and_filter(alignment_file,custom_id=True, extension_length=500, source_fasta = source_fasta)
#alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file,custom_id=True)
# %%
coord_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_main/gene_ref/protein_coding_sequences_sliced.txt'
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/zoonomia447/hg38.phyloP447way.bw'
phylop_447=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', custom_id=True, extension_length=500)
#phylop_447=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', custom_id=True)
#%%
mean_phylop_447=mapper.normalise(alignment=alignment_extended, mapped_data=phylop_447)
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
from ma_mapper import custom_cmap
plots.plot_experimental(heatmap=True,
    data = [phylop_447,], 
    alignment=alignment_extended, 
    #show_alignment=True,
    #alignment_col='dna',
    #show_alignment_colbar=True,
    hm_plot_title =f'phyloP of coding sequences (first/last 500 bp)',
    #hm_plot_title =f'phyloP of coding sequences 50bp from start codons',
    hm_ylabel ='sequences',
    hm_title_fs=12,
    heatmap_yhighlight=[[1000,1000]],
    heatmap_yhighlight_col=['black'],
    heatmap_yhighlight_alpha=[0.5],
    heatmap_color=[custom_cmap.vlag_r_mpl,],
    heatmap_mode='overlay', 
    vlim = [[-1,1],[0,0.0005],], 
    hm_interpolation =  None,
    opacity = 0.8, 
    aggregated=True, 
    aggregated_data=[mean_phylop_447], 
    agg_colset=['grey','grey','red'],
    agg_yhighlight=[[1000,1000]],
    agg_yhighlight_col=['black'],
    agg_yhighlight_alpha=[0.5],
    agg_ylim=[[-0.2,6]],
    #agg_titles=['phyloP_447',], 
    agg_ylabel=['mean phyloP'],
    agg_ylabel_fs=8,
    agg_xlabel='position (bp)',
    colorbar=True,
    colorbar_steps = [0.1,0.0001],  
    #figsize=[100,60],
    agg_major_tick=200,
    agg_h =10,
    #xlim=[450,550],
    #figsize=[60,60]
    )
# %%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
from ma_mapper import custom_cmap
plots.plot_experimental(heatmap=True,
    data = [phylop_447,], 
    alignment=alignment_extended, 
    hm_interpolation =  None,
    #show_alignment=True,
    #alignment_col='dna',
    #show_alignment_colbar=True,
    #hm_plot_title =f'phyloP of coding sequences (first/last 500 bp)',
    hm_plot_title =f'phyloP of coding sequences 50bp from start codons',
    hm_ylabel ='sequences',
    hm_title_fs=12,
    heatmap_yhighlight=[[1000,1000]],
    heatmap_yhighlight_col=['black'],
    heatmap_yhighlight_alpha=[0.5],
    heatmap_color=[custom_cmap.vlag_r_mpl,],
    heatmap_mode='overlay', 
    vlim = [[-1,1],[0,0.0005],], 
    opacity = 0.8, 
    aggregated=True, 
    aggregated_data=[mean_phylop_447], 
    agg_colset=['grey','grey','red'],
    agg_yhighlight=[[1000,1000]],
    agg_yhighlight_col=['black'],
    agg_yhighlight_alpha=[0.5],
    agg_ylim=[[-0.2,6]],
    #agg_titles=['phyloP_447',], 
    agg_ylabel=['mean phyloP'],
    agg_xlabel='position (bp)',
    agg_ylabel_fs=8,
    colorbar=True,
    colorbar_steps = [0.1,0.0001],  
    #figsize=[100,60],
    agg_major_tick=10,
    agg_h =20,
    xlim=[450,550],
    figsize=[60,30]
    )
# %%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
from ma_mapper import custom_cmap
plots.plot_experimental(heatmap=True,
    #data = [phylop_447,], 
    alignment=alignment_extended, 
    hm_plot_title =f'Alignment of coding sequences 50bp from start codons',
    hm_ylabel ='sequences',
    hm_xlabel = 'position (bp)',
    hm_title_fs=12,
    heatmap_yhighlight=[[1000,1000]],
    heatmap_yhighlight_col=['black'],
    heatmap_yhighlight_alpha=[0.5],
    heatmap_color=[custom_cmap.vlag_r_mpl,],
    show_alignment=True,
    alignment_col='dna',
    show_alignment_colbar=True,
    #heatmap_color=[custom_cmap.vlag_r_mpl,],
    heatmap_mode='overlay', 
    vlim = [[-0.5,0.5],[0,0.0005],], 
    opacity = 0.8, 
    colorbar=True,
    agg_major_tick=10,
    xlim=[450,550],
    figsize = [60,30]
    )
# %%
# %%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
from ma_mapper import custom_cmap
plots.plot_experimental(heatmap=True,
    data = [phylop_447,], 
    alignment=alignment_extended, 
    #show_alignment=True,
    #alignment_col='dna',
    #show_alignment_colbar=True,
    hm_plot_title =f'phyloP of coding sequences 50bp from termination codons',
    hm_ylabel ='sequences',
    hm_title_fs=12,
    heatmap_yhighlight=[[1000,1000]],
    heatmap_yhighlight_col=['black'],
    heatmap_yhighlight_alpha=[0.5],
    heatmap_color=[custom_cmap.vlag_r_mpl,],
    heatmap_mode='spread_horizontal', 
    vlim = [[-1,1],[0,0.0005],], 
    hm_interpolation =  None,
    opacity = 0.8, 
    agg_ylabel=['mean phyloP'],
    agg_ylabel_fs=8,
    agg_xlabel='position (bp)',
    aggregated=True, 
    aggregated_data=[mean_phylop_447], 
    agg_colset=['grey','grey','red'],
    agg_ylim=[[-1.2,6]],
    #agg_titles=['phyloP_447',], 
    colorbar=True,
    colorbar_steps = [0.1,0.0001],  
    #figsize=[100,60],
    agg_major_tick=10,
    agg_h =20,
    xlim=[1450,1550],
    figsize=[60,30]
    )
# %%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
from ma_mapper import custom_cmap
plots.plot_experimental(heatmap=True,
    #data = [phylop_447,], 
    alignment=alignment_extended, 
    show_alignment=True,
    alignment_col='dna',
    show_alignment_colbar=True,
    hm_plot_title =f'Alignment of coding sequences 50bp from termination codons',
    hm_ylabel ='sequences',
    hm_xlabel = 'position (bp)',
    hm_title_fs=12,
    #heatmap_color=[custom_cmap.vlag_r_mpl,],
    heatmap_mode='spread_horizontal', 
    vlim = [[-0.5,0.5],[0,0.0005],], 
    opacity = 0.8, 
    colorbar=True,
    agg_major_tick=10,
    xlim=[1450,1550],
    figsize = [60,30]
    )
# %%
