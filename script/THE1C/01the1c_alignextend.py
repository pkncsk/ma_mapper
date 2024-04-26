#%% extract sequence
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
subfamily = ['THE1C']
from ma_mapper import sequence_alignment
coord_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/coord_internal_id/'+subfamily[0]+'.txt'
sequence_alignment.extract_subfamily_coord(subfamily, output_filepath= coord_file)
#%% alignment
import sys
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import sequence_alignment
subfamily = ['THE1C']
source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_fasta/hg38.fa'
metadata = coord_file
fasta_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/alignment/'+subfamily[0]+'.fasta'
sequence_alignment.fetch_sequence(metadata,source_fasta,output_filepath =fasta_file, save_to_file= True,custom_id= True)
sequence_alignment.mafft_align(fasta_file, nthread = 40)
#%% parse alignment
import sys
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
subfamily = ['THE1C']
from ma_mapper import mapper
alignment_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/alignment/'+subfamily[0]+'.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered)
#%%
from ma_mapper import plot
plot.overlay_plot([alignment_filtered])
#%%
plot.overlay_plot([alignment_filtered], metadata= metadata_age)
#%%
te_age = metadata_age.te_age
age_unique=te_age.sort_values().unique()
colors_hex = ['#000000','#252525','#525252',  '#737373',  '#969696',  '#bdbdbd', '#d9d9d9',  '#f0f0f0',  '#ffffff']


#%%
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.colors import BoundaryNorm
import os
import sys
import numpy as np
age_unique = [-0.5, 0.  ,  6.7 ,  9.06, 15.76, 20.19, 29.44, 43.2 , 73.8 , 76.  ]
plt.rcParams['figure.dpi'] = 600
nucleotide_color =['grey','green','yellow','red','blue']
nucleotide_labels = ['gap', 'A', 'C', 'T', 'G']
align_cmap = ListedColormap(nucleotide_color)
ai_blue = ['#000000','#252525','#525252',  '#737373',  '#969696',  '#bdbdbd', '#d9d9d9',  '#f0f0f0',  '#ffffff']
ai_blue_cmap = ListedColormap(ai_blue)
norm = BoundaryNorm(age_unique, len(ai_blue))
fig = plt.figure(figsize=(8,5))
grid = fig.add_gridspec(nrows = 12, ncols = 30, hspace=0)
ax0 = fig.add_subplot(grid[0:10,0:21])
ax1 =fig.add_subplot(grid[0:10,21])
ax2 =fig.add_subplot(grid[5:10,23:24])
ax3 = fig.add_subplot(grid[0:3,23:24])
heatmap_core = ax0.imshow(alignment_filtered, aspect = 'auto',cmap= align_cmap, interpolation='nearest', vmax=5)
ax0.set_yticks([])
ax0.set_title('THE1C alignment')
row_annot=ax1.imshow(te_age.values.reshape(-1, 1), aspect = 'auto',cmap= ai_blue_cmap, interpolation='nearest', norm=norm)
ax1.set_xticks([])
ax1.set_yticks([])
cbar1 = fig.colorbar(row_annot, cax=ax2)
cbar1.set_ticks([(age_unique[i] + age_unique[i+1]) / 2 for i in range(len(age_unique) - 1)])
cbar1.set_ticklabels(age_unique[1:], fontsize = 'small')
cbar1.ax.set_title('te_age', fontsize ='small')
ax2.set_xticks([])
cbar2 = fig.colorbar(heatmap_core, cax=ax3, ticks=[0.5, 1.5, 2.5, 3.5, 4.5], orientation='vertical')
cbar2.set_ticklabels(nucleotide_labels, fontsize='small')
cbar2.ax.set_title('nucleotide', fontsize='small')
plt.tight_layout()
plt.show()
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
subfamily = ['THE1C']
import logging
logging.getLogger('matplotlib').setLevel(logging.WARNING)
#%% from age_div table
species_reference = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/combined_age_div/combined_age_and_div.txt'
main_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
age_div_table = pd.read_csv(species_reference, sep='\t')
subfam_table=age_div_table[age_div_table.repName.isin(subfamily)]
subfam_table = subfam_table[subfam_table.genoName.isin(main_chr)]
subfam_table['id'] = subfam_table.repName+'_'+subfam_table.internal_id.astype(str)
#%%
subfam_age=subfam_table[['id','te_age','te_div']].drop_duplicates()
metadata_with_te_age=metadata_aligned_filtered.merge(subfam_age, on = 'id', how ='left')
#%%
age_colorcode = {0:'#bf6f93',6.7:'#79b743',9.06:'#b262cb',15.76:'#cea240',20.19:'#6566d1',29.44:'#738139',43.2:'#ca489e',73.8:'#51ad7b',76:'#cf4766',82:'#45b0cf',90:'#d25337',96:'#7981c7',105:'#b97446'}
metadata_with_te_age[['subfam', 'subfam_order']] = metadata_with_te_age['id'].str.split('_', expand=True)
row_color_age=metadata_with_te_age.te_age.map(age_colorcode)
map_color = ['grey','green','yellow','red','blue']
custom_cmap = LinearSegmentedColormap.from_list('Custom', map_color, len(map_color))

#%%
framed_alignment=pd.DataFrame(aligned_filtered)
graphical_object=sns.clustermap(framed_alignment,row_colors=row_color_age, row_cluster=False, col_cluster=False, cmap =  custom_cmap,xticklabels =aligned_filtered.shape[1]-1, yticklabels = 500, annot = False)
graphical_object.fig.subplots_adjust(left=0.05)
graphical_object.ax_cbar.set_position((0.15, .08, .02, .4))
graphical_object.ax_cbar.set_ylabel('allele')
graphical_object.ax_cbar.set_yticklabels(['gap','','A','','C','','T','','G'])
col = graphical_object.ax_col_dendrogram.get_position()
graphical_object.ax_col_dendrogram.set_position([col.x0, col.y0, col.width*0.25, col.height*0.25])
for label in metadata_with_te_age.te_age.unique():
    graphical_object.ax_row_dendrogram.bar(0, 0, color=age_colorcode[label], label=label, linewidth=0)
l1 = graphical_object.ax_row_dendrogram.legend(title='te_age', loc="upper right", bbox_to_anchor=(0.2, 0.8), bbox_transform=gcf().transFigure)

graphical_object.ax_heatmap.set_title(subfamily[0]+" alignment")
plt.setp(graphical_object.ax_heatmap.set_xlabel("position (bp)"))
plt.setp(graphical_object.ax_heatmap.set_ylabel("sequences"))
plt.show()
#%% EXTENSION
coord_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/coord_internal_id/'+subfamily[0]+'.txt'
#%%
metadata_df = pd.read_csv(coord_file, sep='\t')
original_order = metadata_df.iloc[:,4].unique()
#%% find border for each id
chrom_list = []
low_border_list = []
high_border_list = []
strand_list = []
for uniq_meta_id in original_order:
    metadata_by_id = metadata_df[metadata_df.id == uniq_meta_id]
    chrom_list.append(metadata_by_id.iloc[:,0].unique()[0])
    low_border_list.append(min(metadata_by_id.iloc[:,1]))
    high_border_list.append(max(metadata_by_id.iloc[:,2]))
    strand_list.append(metadata_by_id.iloc[:,3].unique()[0])
#%% make new metadata 
source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/hg38.fa'
#%%
temp_dict = {'chrom':chrom_list,'start':low_border_list,'end':low_border_list,'strand':strand_list,'id':original_order}
low_border_metadata = pd.DataFrame(temp_dict)
low_border_metadata.start = low_border_metadata.start-500
low_border_records=sequence_alignment.fetch_sequence(low_border_metadata,source_fasta, custom_id= False)

# %%
temp_dict = {'chrom':chrom_list,'start':high_border_list,'end':high_border_list,'strand':strand_list,'id':original_order}
high_border_metadata = pd.DataFrame(temp_dict)
high_border_metadata.end = high_border_metadata.end+500
high_border_records=sequence_alignment.fetch_sequence(high_border_metadata,source_fasta, custom_id= False)
# %%
front_list = []
back_list = []
for idx, strand in enumerate(strand_list):
    if strand == '+':
        front_list.append(low_border_records[idx])
        back_list.append(high_border_records[idx])
    else:
        front_list.append(high_border_records[idx])
        back_list.append(low_border_records[idx])
    
# %%
front_parsed = mapper.parse_alignment(front_list, save_to_file= False)
back_parsed = mapper.parse_alignment(back_list, save_to_file= False)
#%%
filters=mapper.create_filter(aligned_parsed)
row_filter = filters[0]
col_filter = filters[1]
aligned_col_filtered=aligned_parsed[np.ix_(range(len(aligned_parsed)),col_filter)]
#%%
front_parsed_sorted = []
back_parsed_sorted = []
for idx, row in metadata_aligned.iterrows():
    front_parsed_sorted.append(front_parsed[np.where(original_order == row.id)[0][0]])
    back_parsed_sorted.append(back_parsed[np.where(original_order == row.id)[0][0]])
# %%
fused_parsed = []
for i in range(len(aligned_parsed)):
    fused_parsed.append(np.concatenate((front_parsed_sorted[i], aligned_col_filtered[i], back_parsed_sorted[i])))
fused_parsed = np.array(fused_parsed)
#%%
fused_parsed_row_filtered = fused_parsed[np.ix_(row_filter,range(fused_parsed.shape[1]))]
#%%
metadata_aligned_filtered=metadata_aligned.iloc[row_filter,:]
#%% from age_div table
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
metadata_with_te_age[['subfam', 'subfam_order']] = metadata_with_te_age['id'].str.split('_', expand=True)
row_color_age=metadata_with_te_age.te_age.map(age_colorcode)

map_color = ['grey','green','yellow','red','blue']
custom_cmap = LinearSegmentedColormap.from_list('Custom', map_color, len(map_color))

#%%
framed_alignment=pd.DataFrame(fused_parsed_row_filtered)
graphical_object=sns.clustermap(framed_alignment,row_colors=row_color_age, row_cluster=False, col_cluster=False, cmap =  custom_cmap,xticklabels =fused_parsed_row_filtered.shape[1]-1, yticklabels = 500, annot = False)
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

graphical_object.ax_heatmap.set_title(subfamily[0]+" alignment")
plt.setp(graphical_object.ax_heatmap.set_xlabel("position (bp)"))
plt.setp(graphical_object.ax_heatmap.set_ylabel("sequences"))
plt.show()
#%%