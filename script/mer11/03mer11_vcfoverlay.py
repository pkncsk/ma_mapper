#%%
import sys
import pandas as pd
import numpy as np
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import mapper
from ma_mapper import fetch_data
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import gcf
#%%
input_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11.fasta.aligned'
aligned_parsed = mapper.parse_alignment(input_filepath, save_to_file= False)
metadata_aligned = mapper.extract_metadata_from_alignment(input_filepath)
#%%
metadata_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11_coord_with_id.txt'
vcf_mapped=fetch_data.fetch_vcf(metadata_input= metadata_filepath, vcf_input='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/vcf-gnomad/', vcf_format = 'gnomad', custom_id= True)
# %%
metadata_df = pd.read_csv(metadata_filepath, sep='\t')
original_order = metadata_df.iloc[:,4].unique()
#%%
vcf_mapped_sorted = []
for idx, row in metadata_aligned.iterrows():
    vcf_mapped_sorted.append(vcf_mapped[np.where(original_order == row.id)[0][0]])
#%%
filters=mapper.create_filter(aligned_parsed)
# %%
aligned_vcf_overlay=mapper.map_data(vcf_mapped_sorted, aligned_parsed, filters = filters)
#%%
row_filter = filters[0]
col_filter = filters[1]
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
#%%

framed_alignment=pd.DataFrame(aligned_vcf_overlay)
graphical_object=sns.clustermap(framed_alignment,row_colors=row_colors, row_cluster=False, col_cluster=False, cmap =  "Blues",xticklabels =aligned_vcf_overlay.shape[1]-1, yticklabels = 500, annot = False, vmax=0.0005)
graphical_object.fig.subplots_adjust(left=0.05)
graphical_object.ax_cbar.set_position((0.1, .08, .02, .4))
graphical_object.ax_cbar.set_ylabel('normalised_alt_allele_freq')
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
graphical_object.ax_heatmap.set_title("MER11 VCF overlay")
plt.setp(graphical_object.ax_heatmap.set_xlabel("position (bp)"))
plt.setp(graphical_object.ax_heatmap.set_ylabel("sequences"))
plt.show()
# %%
import scipy
#%%
distance_matrix = scipy.spatial.distance.pdist(aligned_vcf_overlay, metric='euclidean')
linkage_array = scipy.cluster.hierarchy.linkage(distance_matrix, method = 'ward')
new_order = scipy.cluster.hierarchy.leaves_list(linkage_array)
aligned_vcf_overlay_sorted = []
for i in new_order:
    aligned_vcf_overlay_sorted.append(aligned_vcf_overlay[i])
aligned_vcf_overlay_sorted = np.array(aligned_vcf_overlay_sorted)
# %%
metadata_with_te_age_sorted=metadata_with_te_age.iloc[new_order]
# %%
age_colorcode = {0:'#bf6f93',6.7:'#79b743',9.06:'#b262cb',15.76:'#cea240',20.19:'#6566d1',29.44:'#738139',43.2:'#ca489e',73.8:'#51ad7b',76:'#cf4766',82:'#45b0cf',90:'#d25337',96:'#7981c7',105:'#b97446'}
subfam_colorcode={'MER11A': 'red', 'MER11B':'blue', 'MER11C':'yellow'}
metadata_with_te_age_sorted[['subfam', 'subfam_order']] = metadata_with_te_age_sorted['id'].str.split('_', expand=True)
row_color_subfam_sorted=metadata_with_te_age_sorted.subfam.map(subfam_colorcode)
row_color_age_sorted=metadata_with_te_age_sorted.te_age.map(age_colorcode)
row_colors_sorted=pd.DataFrame({'subfam':row_color_subfam_sorted,'te_age':row_color_age_sorted}).reset_index(drop=True)
# %%
graphical_object_sorted=sns.clustermap(pd.DataFrame(aligned_vcf_overlay_sorted),row_colors=row_colors_sorted, row_cluster=False, col_cluster=False, cmap =  "Blues",xticklabels =aligned_vcf_overlay_sorted.shape[1]-1, yticklabels = 500, annot = False, vmax=0.0005)
graphical_object_sorted.fig.subplots_adjust(left=0.05)
graphical_object_sorted.ax_cbar.set_position((0.1, .08, .02, .4))
graphical_object_sorted.ax_cbar.set_ylabel('normalised_alt_allele')
col = graphical_object_sorted.ax_col_dendrogram.get_position()
graphical_object_sorted.ax_col_dendrogram.set_position([col.x0, col.y0, col.width*0.25, col.height*0.25])
#graphical_object.ax_row_dendrogram.set_visible(False)
#graphical_object.ax_col_dendrogram.set_visible(False)
for label in metadata_with_te_age_sorted.te_age.unique():
    graphical_object_sorted.ax_row_dendrogram.bar(0, 0, color=age_colorcode[label], label=label, linewidth=0)
l1 = graphical_object_sorted.ax_row_dendrogram.legend(title='te_age', loc="upper right", bbox_to_anchor=(0.2, 0.8), bbox_transform=gcf().transFigure)
for label in metadata_with_te_age_sorted.subfam.unique():
    graphical_object_sorted.ax_col_dendrogram.bar(0, 0, color=subfam_colorcode[label], label=label, linewidth=0)
l2 = graphical_object_sorted.ax_col_dendrogram.legend(title='subfamily', loc="upper right", bbox_to_anchor=(0.205, 0.6), bbox_transform=gcf().transFigure)
graphical_object_sorted.ax_heatmap.set_title("MER11 VCF overlay clustered")
plt.setp(graphical_object_sorted.ax_heatmap.set_xlabel("position (bp)"))
plt.setp(graphical_object_sorted.ax_heatmap.set_ylabel("sequences"))
plt.show()
# %%
