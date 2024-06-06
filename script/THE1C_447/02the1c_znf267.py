#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import extract_bigwig
from ma_mapper import extract_bam
from ma_mapper import mapper
#%%
subfamily = ['THE1C']
alignment_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/alignment/'+subfamily[0]+'.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered)
# %%
coord_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/coord_internal_id/'+subfamily[0]+'.txt'
#%%
bam_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/KZFP-bam_hg38/znf267.sorted.bam'
#%%
bam_forward=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_forward')
bam_reverse=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_reverse')
bam_min=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_min')

# %%
from ma_mapper import plot
plot.overlay_plot([alignment_filtered,bam_forward, bam_reverse], metadata= metadata_age, nucleotide_color='white',data_vmax = [0.1,0.1,0.1],data_cmap=['Blues','Reds'] plot_title='ChIPexo_znf267', show_data_legend=True,data_label=['ChIP_forward','ChIP_reverse'])
#%%
plot.overlay_plot([alignment_filtered,bam_min], metadata= metadata_age, nucleotide_color='white',data_vmax = [0.1,0.1,0.1],data_cmap=['Greens'], plot_title='ChIPexo_znf267', show_data_legend=True,data_label=['ChIP_minimum'])