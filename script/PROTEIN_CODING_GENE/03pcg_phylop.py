#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import extract_bigwig
from ma_mapper import extract_bam
from ma_mapper import mapper
#%%
subfamily = ['pcg_sliced']
alignment_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/alignment/'+subfamily[0]+'.fasta'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered)
# %%
coord_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/coord_internal_id/'+subfamily[0]+'.txt'
bigwig_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_bigwig/241-mammalian-2020v2.bigWig'
#%%
phylop=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig')

# %%
from ma_mapper import plot
plot.overlay_plot([alignment_filtered,phylop], metadata= metadata_age,data_vmin=[-1.0],data_vmax=[1.0], nucleotide_color='white',data_cmap=['RdBu'], plot_title='pcg_sliced:Zoonomia_phyloP', show_data_legend=True,data_label=['phyloP'])
#%%
import matplotlib.pyplot as plt
import numpy as np
normaliser = np.count_nonzero(alignment_filtered, axis=0)
fig, ax = plt.subplots(figsize=(10,3))
ax.fill_between(range(phylop.shape[1]), np.nansum(phylop, axis=0)/phylop.shape[0], color = 'grey')
ax.margins(x=0, y=0)
ax.set_ylim(-0.5,0.5)
# %%
