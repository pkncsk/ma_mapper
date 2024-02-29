#%%
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
metadata_aligned['original_order'] = metadata_aligned.index
#%%
metadata_aligned_a = metadata_aligned[metadata_aligned.id.str.contains('MER11A')]
metadata_aligned_b = metadata_aligned[metadata_aligned.id.str.contains('MER11B')]
metadata_aligned_c = metadata_aligned[metadata_aligned.id.str.contains('MER11C')]
#%% rearrange
metadata_sorted=pd.concat([metadata_aligned_a,metadata_aligned_b,metadata_aligned_c])
aligned_sorted = []
for idx, row in metadata_sorted.iterrows():
    aligned_sorted.append(aligned_parsed[row.original_order])
aligned_sorted=np.array(aligned_sorted)
#%%
filters=mapper.create_filter(aligned_sorted)
# %%
row_filter = filters[0]
col_filter = filters[1]
aligned_sorted_filtered=aligned_sorted[np.ix_(row_filter,col_filter)]

#%%
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce
from matplotlib.colors import LinearSegmentedColormap

#%%
map_color = ['grey','green','yellow','red','blue']
custom_cmap = LinearSegmentedColormap.from_list('Custom', map_color, len(map_color))

metadata_sorted_filtered=metadata_sorted.iloc[row_filter,:]
subfam_colorcode={'MER11A': 'red', 'MER11B':'blue', 'MER11C':'yellow'}
metadata_sorted_filtered[['subfam', 'subfam_order']] = metadata_sorted_filtered['id'].str.split('_', expand=True)
row_color=metadata_sorted_filtered.subfam.map(subfam_colorcode)
graphical_object=sns.clustermap(aligned_sorted_filtered,row_colors=row_color.to_numpy(), row_cluster=False, col_cluster=False, cmap =  custom_cmap,xticklabels =aligned_sorted_filtered.shape[1]-1, yticklabels = 500, annot = False)
graphical_object.cax.set_visible(False)
graphical_object.ax_row_dendrogram.set_visible(False)
graphical_object.ax_col_dendrogram.set_visible(False)
graphical_object.ax_heatmap.set_title("MER11 alignment")
plt.setp(graphical_object.ax_heatmap.set_xlabel("position (bp)"))
plt.setp(graphical_object.ax_heatmap.set_ylabel("sequences"))
plt.show()

# %%
