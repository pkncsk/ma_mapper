#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import extract_bigwig, mapper, sequence_alignment
from ma_mapper import custom_cmap
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_main as config
import compress_pickle
#%%
subfamily = 'THE1C'
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
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
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/zoonomia447/hg38.phyloP447way.bw'
phylop_447=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', custom_id=True)
# %%
mean_phylop=mapper.normalise(alignment=alignment_filtered, mapped_data=phylop_447)
#%%
from ma_mapper import plots
plots.plot_experimental(data = [phylop_447], heatmap_color=[custom_cmap.vlag_r_mpl], vlim = [[-0.5,0.5]], opacity = 0.5, colorbar = True, colorbar_steps = [0.1])
#%%
mean_phylop=mapper.normalise(alignment=alignment_filtered, mapped_data=phylop_447)
import importlib 
importlib.reload(plots)
plots.plot_experimental(data = [phylop_447], heatmap_color=[custom_cmap.vlag_r_mpl], vlim = [[-0.5,0.5]] , opacity = 0.5,colorbar=True, colorbar_steps = [0.1],aggregated = True, agg_colset = ['grey'], aggregated_data=[mean_phylop], agg_ylim=[[-1,1]])
# %%
alignment_extend, metadata_filtered= mapper.parse_and_filter(alignment_file, extension_length=500, source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_fasta/hg38.fa')
phylop_447_ext=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', extension_length=500, source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_fasta/hg38.fa', custom_id=True)
mean_phylop_ext=mapper.normalise(alignment=alignment_extend, mapped_data=phylop_447_ext)
# %%
import importlib 
importlib.reload(plots)
plots.plot_experimental(data = [phylop_447_ext], heatmap_color=[custom_cmap.vlag_r_mpl], vlim = [[-0.5,0.5]] ,opacity = 0.5,colorbar=True, colorbar_steps = [0.1],aggregated = True, agg_colset = ['grey'], aggregated_data=[mean_phylop_ext])
# %%
