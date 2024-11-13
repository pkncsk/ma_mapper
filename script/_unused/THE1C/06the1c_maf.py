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
maf_dir = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/241genomes/241-mammalian-2020v2b.maf'
species_age_table = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/species241_info.tsv'
#%%
# IMPORTANT:
# count_arg is now set to count_arg='common' by default this one counts the common base and then finds the percentage of it against total alignment
# human_ref = counts same the base as human reference then find a percentage
# coverage/total_raw = counts numbers of alignment
# common_raw = finds a raw count of the common base
# if older results need to be reproduced, use:
#maf=mapper.map_and_overlay(alignment_file, coord_file, maf_dir, data_format='maf', separated_maf = True, age_arg='calibrate', species_age_table_file=species_age_table, count_arg='human_ref')  
maf=mapper.map_and_overlay(alignment_file, coord_file, maf_dir, data_format='maf', separated_maf = True, age_arg='calibrate', species_age_table_file=species_age_table)

# %%
from ma_mapper._unused import plot
plot.overlay_plot([alignment_filtered,maf], metadata= metadata_age,data_vmin=[0],data_vmax=[1.0], nucleotide_color='white',data_cmap=['viridis'], plot_title='THE1C:Zoonomia_maf', show_data_legend=True,data_label=['maf'])

# %%
