#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_main as config
from ma_mapper import mapper
from ma_mapper import custom_cmap
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.stats import false_discovery_control
def fdr_control_with_nans(p_values, method='bh'):
    p_values = np.asarray(p_values)
    valid_mask = ~np.isnan(p_values)
    valid_p_values = p_values[valid_mask]
    adj_p_values_valid = false_discovery_control(valid_p_values, method=method)
    adj_p_values = np.full(p_values.shape, np.nan)
    adj_p_values[valid_mask] = adj_p_values_valid
    return adj_p_values
age_ref_table = pd.DataFrame(config.age_ref_table_template)
#%%
subfamily='MER11'
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
#%%
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file,col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
age_table_list = [f'{config.te_age_folder}/MER11A.txt',
             f'{config.te_age_folder}/MER11B.txt',
             f'{config.te_age_folder}/MER11C.txt']
age_df_list = []
for age_tbl in age_table_list:
    age_df_list.append(pd.read_csv(age_tbl, sep='\t'))
age_df=pd.concat(age_df_list)
internal_id_tbl_list = [f'{config.internal_id_folder}/MER11A.txt',
             f'{config.internal_id_folder}/MER11B.txt',
             f'{config.internal_id_folder}/MER11C.txt']
internal_id_df_list = []
for internal_id_tbl in internal_id_tbl_list:
    internal_id_df_list.append(pd.read_csv(internal_id_tbl, sep='\t').sort_values('rmsk_index'))
internal_id_sort=pd.concat(internal_id_df_list)
te_age_internal_id=internal_id_sort.merge(age_df, on='internal_id', how='left')
#%%
age_default_id = pd.DataFrame()
age_default_id['internal_id'] = subfamily + '_' + te_age_internal_id.index.astype(str)
age_default_id['group'] = te_age_internal_id.internal_id
age_default_id['te_age'] = te_age_internal_id['te_age']
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered, age_table=age_default_id)

#%%
###############IMPORTANT####################
#filter NA
metadata_age=metadata_age[~metadata_age['te_age'].isna()]
noNA_indices=metadata_age.index
metadata_age=metadata_age.reset_index()
alignment_filtered=alignment_filtered[noNA_indices]
#%%
import numpy as np
metadata_age['len'] = metadata_age.end.astype(int) - metadata_age.start.astype(int)
age_subgroups = np.unique(metadata_age['te_age'].sort_values())
age_subgroup = {subgroup: num for num, subgroup in enumerate(age_subgroups)}
age_anno=metadata_age['te_age'].map(age_subgroup)
metadata_age['subfam'] = metadata_age['group'].str.split('_').str[0]
subfam_subgroups = np.unique(metadata_age['subfam'].astype(str))
subfam_subgroup = {subgroup: num for num, subgroup in enumerate(subfam_subgroups)}
subfam_anno=metadata_age['subfam'].map(subfam_subgroup)
#%%
coord_file=mapper.extract_metadata_from_alignment(alignment_file)
#%%
bam_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/KZFP-bam_hg38/znf808.sorted.bam'
znf808=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_sum', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
znf808=znf808[noNA_indices]

#%%
bam_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/KZFP-bam_hg38/znf525.sorted.bam'
znf525=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_sum', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
znf525=znf525[noNA_indices]
#%%
bam_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/KZFP-bam_hg38/znf440.sorted.bam'
znf440=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_sum', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
znf440=znf440[noNA_indices]
#%%
bam_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/KZFP-bam_hg38/znf468.sorted.bam'
znf468=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_sum', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
znf468=znf468[noNA_indices]
#%%
bam_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/KZFP-bam_hg38/znf433.sorted.bam'
znf433=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_sum', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
znf433=znf433[noNA_indices]
#%%
bam_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/KZFP-bam_hg38/znf578.sorted.bam'
znf578=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_sum', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
znf578=znf578[noNA_indices]
#%%
metadata_age_sorted=metadata_age.sort_values(['subfam'])
sorted_indices= metadata_age_sorted.index
znf808_sorted=znf808[sorted_indices]
znf525_sorted=znf525[sorted_indices]
znf440_sorted=znf440[sorted_indices]
znf468_sorted=znf468[sorted_indices]
znf433_sorted=znf433[sorted_indices]
znf578_sorted=znf578[sorted_indices]
alignment_sorted = alignment_filtered[sorted_indices]

subfam_anno_sorted = subfam_anno[sorted_indices]
#%%
sigma = 15
g_radius = None
znf808_blurred = np.array([gaussian_filter1d(row, sigma, radius=g_radius) for row in znf808_sorted])
znf525_blurred = np.array([gaussian_filter1d(row, sigma, radius=g_radius) for row in znf525_sorted])
znf440_blurred = np.array([gaussian_filter1d(row, sigma, radius=g_radius) for row in znf440_sorted])
znf468_blurred = np.array([gaussian_filter1d(row, sigma, radius=g_radius) for row in znf468_sorted])
znf433_blurred = np.array([gaussian_filter1d(row, sigma, radius=g_radius) for row in znf433_sorted])
znf578_blurred = np.array([gaussian_filter1d(row, sigma, radius=g_radius) for row in znf578_sorted])
#mean_min=np.minimum(mean_forward,mean_reverse)
#%%
# %% 
main_minmax = [0,0.20]
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(heatmap=True,
    #data = [blurred_matrix],
    data = [znf808_blurred,znf525_blurred,znf440_blurred,znf468_blurred,znf433_blurred,znf578_blurred], 
    alignment=alignment_filtered, 
    heatmap_color=['Blues','Greens','Purples',custom_cmap.LimeGreens_mpl,'Oranges','Reds'],
    heatmap_mode='overlay', 
    vlim = [main_minmax,main_minmax,main_minmax,main_minmax,main_minmax,main_minmax], 
    opacity = 0.9,
    show_alignment=True, 
    hm_transparency_mode = 'gradient',
        anno_ylabel = 'sequences',
    anno_ylabel_fs=10,
    annotation=True, 
    anno_col = [['red', 'yellow', 'blue']], 
    #anno_title=['TEA-TIME','subfamily'],
    annotation_data=[subfam_anno_sorted],
    anno_cbar_label=[['MER11A','MER11B','MER11C']],
    anno_cbar_title=['subfamily'],  
    anno_cbar_even_pos=-0.19,
    #hm_interpolation = None,
    #aggregated=True, 
    #aggregated_data=[mean_vcf, mean_phylop_447], 
    #agg_colset=['grey','grey','red'],
    #agg_ylim=[[0,0.003],[-1.2,1.2]],
    #agg_ylabel=['average allele freq',None],
    #agg_titles=['gnomAD','phyloP_447',], 
    #colorbar=True,
    #colorbar_steps = [1e-6,],
    #agg_yscale=['log', 'linear'],
    agg_major_tick=100,
    #figsize=[60,20],
    #agg_h=40,   
    hm_plot_title =f'KZFP binding signals on MER11 MSA',
    hm_xlabel = 'position (bp)',
    hm_ylabel = 'sequences',
    )
#%%
bam_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/geo/HNF4A_1.bam'
hnf4=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_sum', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
hnf4=hnf4[noNA_indices]
#%%
bam_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/geo/DUX4_1.bam'
dux4=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_sum', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
dux4=dux4[noNA_indices]
#%%
bam_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/geo/SMAD2.bam'
smad2=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_sum', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
smad2=smad2[noNA_indices]
#%%
bam_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/geo/ZBTB17.bam'
zbtb17=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_sum', custom_id=True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
zbtb17=zbtb17[noNA_indices]
#%%
#%%
metadata_age_sorted=metadata_age.sort_values(['subfam'])
sorted_indices= metadata_age_sorted.index
hnf4_sorted=hnf4[sorted_indices]
dux4_sorted=dux4[sorted_indices]
smad2_sorted=smad2[sorted_indices]
zbtb17_sorted=zbtb17[sorted_indices]
alignment_sorted = alignment_filtered[sorted_indices]

subfam_anno_sorted = subfam_anno[sorted_indices]
#%%
sigma = 15
g_radius = None
hnf4_blurred = np.array([gaussian_filter1d(row, sigma, radius=g_radius) for row in hnf4_sorted])
dux4_blurred = np.array([gaussian_filter1d(row, sigma, radius=g_radius) for row in dux4_sorted])
smad2_blurred = np.array([gaussian_filter1d(row, sigma, radius=g_radius) for row in smad2_sorted])
zbtb17_blurred = np.array([gaussian_filter1d(row, sigma, radius=g_radius) for row in zbtb17_sorted])
#%%
main_minmax = [0,0.20]
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(heatmap=True,
    #data = [blurred_matrix],
    data = [hnf4_blurred,dux4_blurred,smad2_blurred,zbtb17_blurred], 
    alignment=alignment_filtered, 
    heatmap_color=['Purples','Reds','Purples',custom_cmap.LimeGreens_mpl,'Oranges','Reds'],
    heatmap_mode='overlay', 
    vlim = [main_minmax,main_minmax,main_minmax,main_minmax,main_minmax,main_minmax], 
    opacity = 0.9,
    show_alignment=True, 
    hm_transparency_mode = 'gradient',
        anno_ylabel = 'sequences',
    anno_ylabel_fs=10,
    annotation=True, 
    anno_col = [['red', 'yellow', 'blue']], 
    #anno_title=['TEA-TIME','subfamily'],
    annotation_data=[subfam_anno_sorted],
    anno_cbar_label=[['MER11A','MER11B','MER11C']],
    anno_cbar_title=['subfamily'],  
    anno_cbar_even_pos=-0.19,
    #hm_interpolation = None,
    #aggregated=True, 
    #aggregated_data=[mean_vcf, mean_phylop_447], 
    #agg_colset=['grey','grey','red'],
    #agg_ylim=[[0,0.003],[-1.2,1.2]],
    #agg_ylabel=['average allele freq',None],
    #agg_titles=['gnomAD','phyloP_447',], 
    #colorbar=True,
    #colorbar_steps = [1e-6,],
    #agg_yscale=['log', 'linear'],
    agg_major_tick=100,
    #figsize=[60,20],
    #agg_h=40,   
    hm_plot_title =f'TF binding signals on MER11 MSA',
    hm_xlabel = 'position (bp)',
    hm_ylabel = 'sequences',
    )
#%%
n_colors =256
#
#  Creating a gradient from white to lime green over 256 steps
white = np.array([1, 1, 1, 1])  # RGBA for white
limegreen = np.array([0.196, 0.804, 0.196, 1])  # RGBA for limegreen (approximation in normalized [0,1] scale)

# Interpolating between white and limegreen to get 256 colors
rgba_gradient_list = [((1 - t) * white + t * limegreen).tolist() for t in np.linspace(0, 1, n_colors)]

# Adjusting alpha from 0 (transparent) to 1 (opaque) across the gradient
rgba_gradient_alpha_adjusted = [(r, g, b, alpha) for (r, g, b, _), alpha in zip(rgba_gradient_list, np.linspace(0, 1, n_colors))]

# Saving the new list with transparency gradient from 0 to 1
rgba_gradient_alpha_adjusted_path = "/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11/revision/test.txt"
with open(rgba_gradient_alpha_adjusted_path, "w") as file:
    file.write(f"[{', '.join(map(str, rgba_gradient_alpha_adjusted))}]")

rgba_gradient_alpha_adjusted_path
# %%
