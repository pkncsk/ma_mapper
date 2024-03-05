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
subfam_coord_age = subfam_table[['genoName','genoStart','genoEnd','strand']]
subfam_coord_age['id'] = subfam_table.repName+'_'+subfam_table.internal_id.astype(str)
subfam_coord_age['te_age'] = subfam_table.te_age
# %%
subfam_coord_age.to_csv('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11_coord_with_id_age.txt', sep='\t', index= False)
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
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import gcf
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import mapper
#%%
input_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11.fasta.aligned'
aligned_parsed = mapper.parse_alignment(input_filepath, save_to_file= False)
metadata_aligned = mapper.extract_metadata_from_alignment(input_filepath)
metadata_aligned['original_order'] = metadata_aligned.index
#%%
filters=mapper.create_filter(aligned_parsed)
row_filter = filters[0]
col_filter = filters[1]
aligned_filtered=aligned_parsed[np.ix_(row_filter,col_filter)]
metadata_aligned_filtered=metadata_aligned.iloc[row_filter,:]

#%% from age_div table
subfamily = ['MER11A','MER11B','MER11C']
input_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/combined_age_div/combined_age_and_div.txt'
main_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
age_div_table = pd.read_csv(input_filepath, sep='\t')
subfam_table=age_div_table[age_div_table.repName.isin(subfamily)]
subfam_table = subfam_table[subfam_table.genoName.isin(main_chr)]
subfam_table['id'] = subfam_table.repName+'_'+subfam_table.internal_id.astype(str)
#%%
subfam_age=subfam_table[['id','te_age','te_div']].drop_duplicates()
metadata_with_te_age=metadata_aligned_filtered.merge(subfam_age, on = 'id', how ='left')
#%%
age_colorcode = {0:'#bf6f93',6.7:'#79b743',9.06:'#b262cb',15.76:'#cea240',20.19:'#6566d1',29.44:'#738139',43.2:'#ca489e',73.8:'#51ad7b',76:'#cf4766',82:'#45b0cf',90:'#d25337',96:'#7981c7',105:'#b97446'}
subfam_colorcode={'MER11A': 'red', 'MER11B':'blue', 'MER11C':'yellow'}
metadata_with_te_age[['subfam', 'subfam_order']] = metadata_with_te_age['id'].str.split('_', expand=True)
row_color_subfam=metadata_with_te_age.subfam.map(subfam_colorcode)
row_color_age=metadata_with_te_age.te_age.map(age_colorcode)
row_colors=pd.DataFrame({'subfam':row_color_subfam,'te_age':row_color_age})

map_color = ['grey','green','yellow','red','blue']
custom_cmap = LinearSegmentedColormap.from_list('Custom', map_color, len(map_color))

#%%
framed_alignment=pd.DataFrame(aligned_filtered)
graphical_object=sns.clustermap(framed_alignment,row_colors=row_colors, row_cluster=False, col_cluster=False, cmap =  custom_cmap,xticklabels =aligned_filtered.shape[1]-1, yticklabels = 500, annot = False)
graphical_object.fig.subplots_adjust(left=0.05)
graphical_object.ax_cbar.set_position((0.15, .08, .02, .4))
graphical_object.ax_cbar.set_ylabel('allele')
graphical_object.ax_cbar.set_yticklabels(['gap','','A','','C','','T','','G'])
col = graphical_object.ax_col_dendrogram.get_position()
graphical_object.ax_col_dendrogram.set_position([col.x0, col.y0, col.width*0.25, col.height*0.25])
#graphical_object.ax_row_dendrogram.set_visible(False)
#graphical_object.ax_col_dendrogram.set_visible(False)
for label in metadata_with_te_age.te_age.unique():
    graphical_object.ax_row_dendrogram.bar(0, 0, color=age_colorcode[label], label=label, linewidth=0)
l1 = graphical_object.ax_row_dendrogram.legend(title='te_age', loc="upper right", bbox_to_anchor=(0.2, 0.8), bbox_transform=gcf().transFigure)
for label in metadata_with_te_age.subfam.unique():
    graphical_object.ax_col_dendrogram.bar(0, 0, color=subfam_colorcode[label], label=label, linewidth=0)
l2 = graphical_object.ax_col_dendrogram.legend(title='subfamily', loc="upper right", bbox_to_anchor=(0.205, 0.6), bbox_transform=gcf().transFigure)
graphical_object.ax_heatmap.set_title("MER11 alignment")
plt.setp(graphical_object.ax_heatmap.set_xlabel("position (bp)"))
plt.setp(graphical_object.ax_heatmap.set_ylabel("sequences"))
plt.show()
#%%
metadata_aligned_a = metadata_with_te_age[metadata_with_te_age.id.str.contains('MER11A')]
metadata_aligned_b = metadata_with_te_age[metadata_with_te_age.id.str.contains('MER11B')]
metadata_aligned_c = metadata_with_te_age[metadata_with_te_age.id.str.contains('MER11C')]
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
age_colorcode = {0:'#bf6f93',6.7:'#79b743',9.06:'#b262cb',15.76:'#cea240',20.19:'#6566d1',29.44:'#738139',43.2:'#ca489e',73.8:'#51ad7b',76:'#cf4766',82:'#45b0cf',90:'#d25337',96:'#7981c7',105:'#b97446'}
map_color = ['grey','green','yellow','red','blue']
subfam_colorcode={'MER11A': 'red', 'MER11B':'blue', 'MER11C':'yellow'}
metadata_sorted[['subfam', 'subfam_order']] = metadata_sorted['id'].str.split('_', expand=True)
row_color_subfam=metadata_sorted.subfam.map(subfam_colorcode)
row_color_age=metadata_sorted.te_age.map(age_colorcode)
row_colors=pd.DataFrame({'subfam':row_color_subfam,'te_age':row_color_age}).reset_index(drop=True)
custom_cmap = LinearSegmentedColormap.from_list('Custom', map_color, len(map_color))
#%%
framed_alignment=pd.DataFrame(aligned_sorted_filtered)
graphical_object=sns.clustermap(framed_alignment,row_colors=row_colors, row_cluster=False, col_cluster=False, cmap =  custom_cmap,xticklabels =aligned_filtered.shape[1]-1, yticklabels = 500, annot = False)
graphical_object.fig.subplots_adjust(left=0.05)
graphical_object.ax_cbar.set_position((0.15, .08, .02, .4))
graphical_object.ax_cbar.set_ylabel('allele')
graphical_object.ax_cbar.set_yticklabels(['gap','','A','','C','','T','','G'])
col = graphical_object.ax_col_dendrogram.get_position()
graphical_object.ax_col_dendrogram.set_position([col.x0, col.y0, col.width*0.25, col.height*0.25])
#graphical_object.ax_row_dendrogram.set_visible(False)
#graphical_object.ax_col_dendrogram.set_visible(False)
for label in metadata_with_te_age.te_age.unique():
    graphical_object.ax_row_dendrogram.bar(0, 0, color=age_colorcode[label], label=label, linewidth=0)
l1 = graphical_object.ax_row_dendrogram.legend(title='te_age', loc="upper right", bbox_to_anchor=(0.2, 0.8), bbox_transform=gcf().transFigure)
for label in metadata_with_te_age.subfam.unique():
    graphical_object.ax_col_dendrogram.bar(0, 0, color=subfam_colorcode[label], label=label, linewidth=0)
l2 = graphical_object.ax_col_dendrogram.legend(title='subfamily', loc="upper right", bbox_to_anchor=(0.205, 0.6), bbox_transform=gcf().transFigure)
graphical_object.ax_heatmap.set_title("MER11 alignment")
plt.setp(graphical_object.ax_heatmap.set_xlabel("position (bp)"))
plt.setp(graphical_object.ax_heatmap.set_ylabel("sequences"))
plt.show()
# %%
