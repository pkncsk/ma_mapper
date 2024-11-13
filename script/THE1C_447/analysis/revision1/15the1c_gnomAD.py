#%%
import sys
from turtle import width
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_main as config
from ma_mapper import mapper
from ma_mapper import custom_cmap
import numpy as np
#%%
subfamily = 'THE1C'
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file, col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10,extension_length=500,source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_fasta/hg38.fa')
#alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
#alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file, col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
age_table = f'{config.te_age_folder}/{subfamily}.txt'
#%%
age_df = pd.read_csv(age_table, sep='\t')
internal_id_tbl = f'{config.internal_id_folder}/{subfamily}.txt'
internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
internal_id_sort = internal_id_df.sort_values('rmsk_index')
te_age_internal_id=internal_id_sort.merge(age_df, on='internal_id', how='left')
#%%
age_default_id = pd.DataFrame()
age_default_id['internal_id'] = subfamily + '_' + te_age_internal_id.index.astype(str)
age_default_id['te_age'] = te_age_internal_id['te_age']
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered, age_table=age_default_id)
#%%
coord_file=mapper.extract_metadata_from_alignment(alignment_file)
#%%
vcf_dir = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/vcf-gnomad'

vcf=mapper.map_and_overlay(alignment_file, coord_file, vcf_dir, data_format='vcf', vcf_format='gnomad', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10, extension_length=500)

# %%
#%%
mean_vcf=mapper.normalise(alignment=alignment_filtered, mapped_data=vcf)
#%%
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/zoonomia447/hg38.phyloP447way.bw'
phylop_447=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10,extension_length=500)
mean_phylop_447=mapper.normalise(alignment=alignment_filtered, mapped_data=phylop_447)

#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(heatmap=False,
    show_alignment=False, 
    aggregated=True, 
    aggregated_data=[phylop_447[0],phylop_447[1],phylop_447[4500],phylop_447[9000],phylop_447[9700]], 
    agg_colset=['grey','grey','grey','grey','grey'],
    agg_ylim=[[-3,3],[-3,3],[-3,3],[-3,3],[-3,3]],
    agg_ylabel=[None,None,'phyloP score',None,None],
    agg_ylabel_fs=10,
    agg_xlabel='position (bp)',
    agg_plot_title=['phyloP of individual THE1C TE',None,None,None,None],
    agg_titles=['THE1C_0','THE1C_1','THE1C_4500','THE1C_9000','THE1C_9700'], 
    agg_titles_pos=[1.15,0.5],
    colorbar=True,
    colorbar_steps = [1e-6,],
    #agg_yscale=['log', 'linear'],
    agg_major_tick=50,
    figsize=[50,50],
    #agg_h=40,   
    hm_plot_title =f'phyloP of individual TE',
    )
# %% 
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(heatmap=True,
    data = [vcf], 
    alignment=alignment_filtered, 
    heatmap_color=['Blues',],
    heatmap_mode='overlay', 
    vlim = [[1e-6,1e-5],], 
    opacity = 1,
    hm_interpolation=None,
    show_alignment=False, 
    #aggregated=True, 
    #aggregated_data=[mean_vcf, mean_phylop_447], 
    #agg_colset=['grey','grey','red'],
    #agg_ylim=[[0,0.003],[-1.2,1.2]],
    #agg_ylabel=['average allele freq',None],
    #agg_titles=['gnomAD','phyloP_447',], 
    colorbar=True,
    colorbar_steps = [1e-6,],
    #agg_yscale=['log', 'linear'],
    agg_major_tick=100,
    figsize=[60,20],
    #agg_h=40,   
    hm_plot_title =f'Alternate allele frequency from gnomAD on THE1C MSA',
    hm_xlabel = 'position (bp)',
    hm_ylabel = 'sequences',
    )

#%%
test_set = vcf[vcf > 0]
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,1)

# We can set the number of bins with the *bins* keyword argument.
ax.hist(np.log10(test_set), bins=100)

plt.show()
#%%
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(heatmap=True,
    data = [phylop_447], 
    alignment=alignment_filtered, 
    heatmap_color=[custom_cmap.vlag_r_mpl,],
    hm_interpolation=None,
    heatmap_mode='overlay', 
    vlim = [[-1,1]], 
    opacity = 0.8, 
    hm_plot_title =f'phyloP of THE1C MSA with 500bp flanks',
    hm_ylabel ='sequences',
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    hm_title_fs=12,
    aggregated=True, 
    aggregated_data=[mean_phylop_447], 
    agg_colset=['grey','grey','red'],
    agg_ylim=[[-1.2,1.2]],
    agg_ylabel=['mean phyloP',None],
    agg_ylabel_fs=10,
    #agg_titles=[], 
    colorbar=True,
    colorbar_steps = [0.125,],
    agg_yscale=['log', 'linear'],
    agg_major_tick=200,
    #figsize=[60,20],
    agg_h=10,   
    )
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(heatmap=True,
    data = [phylop_447], 
    alignment=alignment_filtered, 
    heatmap_color=[custom_cmap.vlag_r_mpl,],
    hm_interpolation=None,
    heatmap_mode='overlay', 
    vlim = [[-1,1]], 
    opacity = 0.8, 
    hm_plot_title =f'phyloP of THE1C MSA with 500bp flanks',
    hm_ylabel ='sequences',
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    hm_title_fs=12,
    aggregated=False, 
    aggregated_data=[mean_phylop_447], 
    agg_colset=['grey','grey','red'],
    agg_ylim=[[-1.2,1.2]],
    agg_ylabel=['mean phyloP',None],
    agg_ylabel_fs=10,
    #agg_titles=[], 
    colorbar=True,
    colorbar_steps = [0.125,],
    agg_yscale=['log', 'linear'],
    agg_major_tick=100,
    figsize=[60,20],
    agg_h=15,   
    )
# %%
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(heatmap=True,
    #data = [phylop_447], 
    alignment=alignment_filtered, 
    alignment_col='dna',
    show_alignment=True,
    show_alignment_colbar=True,
    hm_plot_title =f'THE1C MSA with 500bp flanks',
    hm_ylabel ='sequences',
    hm_xlabel = 'position (bp)',
    hm_title_fs=10,
    #heatmap_color=[custom_cmap.vlag_r_mpl,],
    heatmap_mode='overlay', 
    vlim = [[-0.5,0.5]], 
    opacity = 0.5, 
    #aggregated=True, 
    #aggregated_data=[mean_vcf, mean_phylop_447], 
    #agg_colset=['grey','grey','red'],
    #agg_ylim=[[0,0.003],[-1.2,1.2]],
    #agg_ylabel=['average allele freq',None],
    #agg_titles=['gnomAD','phyloP_447',], 
    colorbar=True,
    colorbar_steps = [0.125,],
    agg_yscale=['log', 'linear'],
    agg_major_tick=100,
    figsize=[60,20],
    agg_h=40,   
    )
# %%
subfamily='THE1C'
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file, col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
age_table = f'{config.te_age_folder}/{subfamily}.txt'
age_df = pd.read_csv(age_table, sep='\t')
internal_id_tbl = f'{config.internal_id_folder}/{subfamily}.txt'
internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
internal_id_sort = internal_id_df.sort_values('rmsk_index')
te_age_internal_id=internal_id_sort.merge(age_df, on='internal_id', how='left')
age_default_id = pd.DataFrame()
age_default_id['internal_id'] = subfamily + '_' + te_age_internal_id.index.astype(str)
age_default_id['te_age'] = te_age_internal_id['te_age']
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered, age_table=age_default_id)
coord_file=mapper.extract_metadata_from_alignment(alignment_file)
vcf_dir = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/vcf-gnomad'
vcf=mapper.map_and_overlay(alignment_file, coord_file, vcf_dir, data_format='vcf', vcf_format='gnomad', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
#%%
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/zoonomia447/hg38.phyloP447way.bw'
phylop_447=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
mean_phylop_447=mapper.normalise(alignment=alignment_filtered, mapped_data=phylop_447)
mean_vcf=mapper.normalise(alignment=alignment_filtered, mapped_data=vcf)
#%%
# %% 
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(heatmap=True,
    data = [vcf], 
    alignment=alignment_filtered, 
    heatmap_color=['Blues',],
    heatmap_mode='overlay', 
    vlim = [[1e-6,1e-5],], 
    opacity = 1,
    show_alignment=False, 
    hm_interpolation = None,
    #aggregated=True, 
    #aggregated_data=[mean_vcf, mean_phylop_447], 
    #agg_colset=['grey','grey','red'],
    #agg_ylim=[[0,0.003],[-1.2,1.2]],
    #agg_ylabel=['average allele freq',None],
    #agg_titles=['gnomAD','phyloP_447',], 
    colorbar=True,
    colorbar_steps = [1e-6,],
    #agg_yscale=['log', 'linear'],
    agg_major_tick=50,
    #figsize=[60,20],
    #agg_h=40,   
    hm_plot_title =f'Alternate allele frequency\nfrom gnomAD on THE1C MSA',
    hm_xlabel = 'position (bp)',
    hm_ylabel = 'sequences',
    )
# %% 
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(heatmap=True,
    data = [phylop_447], 
    alignment=alignment_filtered, 
    heatmap_color=['RdBu',],
    heatmap_mode='overlay', 
    vlim = [[-1.0,1.0],], 
    opacity = 1,
    show_alignment=False, 
    hm_interpolation = None,
    #aggregated=True, 
    #aggregated_data=[mean_vcf, mean_phylop_447], 
    #agg_colset=['grey','grey','red'],
    #agg_ylim=[[0,0.003],[-1.2,1.2]],
    #agg_ylabel=['average allele freq',None],
    #agg_titles=['gnomAD','phyloP_447',], 
    colorbar=True,
    colorbar_steps = [0.1,],
    #agg_yscale=['log', 'linear'],
    agg_major_tick=50,
    #figsize=[60,20],
    #agg_h=40,   
    hm_plot_title =f'phyloP from Zoonomia dataset on THE1C MSA',
    hm_xlabel = 'position (bp)',
    hm_ylabel = 'sequences',
    )
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(heatmap=True,
    #data = [phylop_447], 
    alignment=alignment_filtered, 
    alignment_col='dna',
    show_alignment=True,
    show_alignment_colbar=True,
    hm_plot_title =f'THE1C MSA',
    hm_ylabel ='sequences',
    hm_xlabel = 'position (bp)',
    hm_title_fs=10,
    #heatmap_color=[custom_cmap.vlag_r_mpl,],
    heatmap_mode='overlay', 
    vlim = [[-0.5,0.5]], 
    hm_interpolation = None,
    opacity = 0.1, 
    #aggregated=True, 
    #aggregated_data=[mean_vcf, mean_phylop_447], 
    #agg_colset=['grey','grey','red'],
    #agg_ylim=[[0,0.003],[-1.2,1.2]],
    #agg_ylabel=['average allele freq',None],
    #agg_titles=['gnomAD','phyloP_447',], 
    colorbar=True,
    colorbar_steps = [0.125,],
    agg_yscale=['log', 'linear'],
    agg_major_tick=50,
    #figsize=[60,20],
    #agg_h=40,   
    )
# %%
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8,8))
flierprops = dict(marker='.', markersize=3, linestyle='none', markeredgecolor='black', alpha=0.6)

ax.scatter(mean_phylop_447, mean_vcf, alpha=0.5, color='black', marker='.')

from scipy import stats
# First regression line (0 to 43.2)
res = stats.linregress(mean_phylop_447, mean_vcf)
ax.plot(mean_phylop_447, res.intercept + res.slope * mean_phylop_447, color='red', label='Regression Line (0-96)', linewidth=2, alpha=0.5)
ax.text(0.99, 0.95, f'n={len(mean_vcf)}', 
        transform=ax.transAxes, 
        fontsize=12, 
        verticalalignment='bottom', 
        horizontalalignment='right')
ax.text(0.99, 0.89, f'regression line: y={res.slope:.2e}x+{res.intercept:.2e}\nR-square={res.rvalue:.2f}', 
        transform=ax.transAxes, 
        fontsize=12,
        color = 'red', 
        verticalalignment='bottom', 
        horizontalalignment='right')
ax.axvline(x=0, color='grey', linewidth=1, alpha=0.5)
#ax.axhline(y=3.342447251181431, color='blue',linestyle='--', linewidth=1, alpha=0.5)
ax.set_ylim(0,0.003)
ax.set_xlabel('mean phyloP per base', fontsize=12)
ax.set_ylabel('mean alternate allele freqeuncy per base', fontsize=12)
ax.set_title('Mean phyloP vs mean alt allele frequency per base on THE1C MSA', fontsize=14)
#plt.xticks(rotation=45, ha='right')
plt.show()
# %%
