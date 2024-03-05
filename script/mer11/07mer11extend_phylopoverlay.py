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
# %%
metadata_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11_coord_with_id.txt'
#%%
metadata_df = pd.read_csv(metadata_filepath, sep='\t')
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
temp_dict = {'chrom':chrom_list,'start':low_border_list,'end':low_border_list,'strand':strand_list,'id':original_order}
low_border_metadata = pd.DataFrame(temp_dict)
low_border_metadata.start = low_border_metadata.start-500
# %%
temp_dict = {'chrom':chrom_list,'start':high_border_list,'end':high_border_list,'strand':strand_list,'id':original_order}
high_border_metadata = pd.DataFrame(temp_dict)
high_border_metadata.end = high_border_metadata.end+500
#%%
metadata_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11_coord_with_id.txt'
bigwig_mapped=fetch_data.fetch_bigwig(metadata_input= metadata_filepath, bigwig_input='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_bigwig/241-mammalian-2020v2.bigWig', custom_id= True)
#%%
metadata_df = pd.read_csv(metadata_filepath, sep='\t')
original_order = metadata_df.iloc[:,4].unique()
#%%
bigwig_mapped_sorted = []
for idx, row in metadata_aligned.iterrows():
    bigwig_mapped_sorted.append(bigwig_mapped[np.where(original_order == row.id)[0][0]])
#%%
filters=mapper.create_filter(aligned_parsed)
row_filter = filters[0]
col_filter = filters[1]
# %%
aligned_filtered=aligned_parsed[np.ix_(row_filter,col_filter)]
normalisation_mask = np.count_nonzero(aligned_filtered, axis=0)
aligned_bigwig_overlay=mapper.map_data(bigwig_mapped_sorted, aligned_parsed, filters = filters)
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

# %%
low_border_bigwig_mapped=fetch_data.fetch_bigwig(metadata_input= low_border_metadata, bigwig_input='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_bigwig/241-mammalian-2020v2.bigWig', custom_id= True)
# %%
high_border_bigwig_mapped=fetch_data.fetch_bigwig(metadata_input= high_border_metadata, bigwig_input='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_bigwig/241-mammalian-2020v2.bigWig', custom_id= True)
#%%
bigwig_front_list = []
bigwig_back_list = []
for idx, strand in enumerate(strand_list):
    if strand == '+':
        bigwig_front_list.append(low_border_bigwig_mapped[idx])
        bigwig_back_list.append(high_border_bigwig_mapped[idx])
    else:
        bigwig_front_list.append(high_border_bigwig_mapped[idx])
        bigwig_back_list.append(low_border_bigwig_mapped[idx])
#%%
bigwig_front_sorted = []
bigwig_back_sorted = []
for idx, row in metadata_aligned.iterrows():
    if row_filter[idx]:
        bigwig_front_sorted.append(bigwig_front_list[np.where(original_order == row.id)[0][0]])
        bigwig_back_sorted.append(bigwig_back_list[np.where(original_order == row.id)[0][0]])
#%%
fused_bigwig_mapped = []
for i in range(len(aligned_bigwig_overlay)):
    fused_bigwig_mapped.append(np.concatenate((bigwig_front_sorted[i], aligned_bigwig_overlay[i], bigwig_back_sorted[i])))
fused_bigwig_mapped = np.array(fused_bigwig_mapped)
#%%
age_colorcode = {0:'#bf6f93',6.7:'#79b743',9.06:'#b262cb',15.76:'#cea240',20.19:'#6566d1',29.44:'#738139',43.2:'#ca489e',73.8:'#51ad7b',76:'#cf4766',82:'#45b0cf',90:'#d25337',96:'#7981c7',105:'#b97446'}
subfam_colorcode={'MER11A': 'red', 'MER11B':'blue', 'MER11C':'yellow'}
metadata_with_te_age[['subfam', 'subfam_order']] = metadata_with_te_age['id'].str.split('_', expand=True)
row_color_subfam_sorted=metadata_with_te_age.subfam.map(subfam_colorcode)
row_color_age_sorted=metadata_with_te_age.te_age.map(age_colorcode)
row_colors_sorted=pd.DataFrame({'subfam':row_color_subfam_sorted,'te_age':row_color_age_sorted}).reset_index(drop=True)
#%%
graphical_object_sorted=sns.clustermap(pd.DataFrame(fused_bigwig_mapped),row_colors=row_colors_sorted, row_cluster=False, col_cluster=False, cmap =  "vlag_r",xticklabels =500, yticklabels = 500, annot = False, vmin = -1, vmax = 1)
graphical_object_sorted.fig.subplots_adjust(left=0.05)
graphical_object_sorted.ax_cbar.set_position((0.1, .08, .02, .4))
graphical_object_sorted.ax_cbar.set_ylabel('phyloP')
col = graphical_object_sorted.ax_col_dendrogram.get_position()
graphical_object_sorted.ax_col_dendrogram.set_position([col.x0, col.y0, col.width*0.25, col.height*0.25])
#graphical_object.ax_row_dendrogram.set_visible(False)
#graphical_object.ax_col_dendrogram.set_visible(False)
for label in metadata_with_te_age.te_age.unique():
    graphical_object_sorted.ax_row_dendrogram.bar(0, 0, color=age_colorcode[label], label=label, linewidth=0)
l1 = graphical_object_sorted.ax_row_dendrogram.legend(title='te_age', loc="upper right", bbox_to_anchor=(0.2, 0.8), bbox_transform=gcf().transFigure)
for label in metadata_with_te_age.subfam.unique():
    graphical_object_sorted.ax_col_dendrogram.bar(0, 0, color=subfam_colorcode[label], label=label, linewidth=0)
l2 = graphical_object_sorted.ax_col_dendrogram.legend(title='subfamily', loc="upper right", bbox_to_anchor=(0.205, 0.6), bbox_transform=gcf().transFigure)
graphical_object_sorted.ax_heatmap.set_title("MER11 phyloP overlay")
plt.setp(graphical_object_sorted.ax_heatmap.set_xlabel("position (bp)"))
plt.setp(graphical_object_sorted.ax_heatmap.set_ylabel("sequences"))
plt.show()
# %%
border_count = np.empty(500, dtype=np.int)
border_count.fill(len(bigwig_front_list))
normaliser=np.concatenate((border_count,normalisation_mask,border_count))
#%%
fused_bigwig_treated=np.nan_to_num(fused_bigwig_mapped)
sum_fused_bigwig_mapped = np.sum(fused_bigwig_treated, axis = 0)
fused_bigwig_averaged = sum_fused_bigwig_mapped/normaliser
# %%
plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600
fig, ax = plt.subplots(figsize=(10,3))
ax.fill_between(range(len(fused_bigwig_averaged)), fused_bigwig_averaged, color='grey')
ax.margins(x=0, y=0)
ax.set_ylim(-1,1)
ax.set_xlabel('position (bp)')
ax.set_ylabel('normalised phyloP')
ax.set_title('MER11 phyloP overlay')
plt.show()

#%%
import scipy
#%%
aligned_bigwig_treated=np.nan_to_num(aligned_bigwig_overlay)
distance_matrix = scipy.spatial.distance.pdist(aligned_bigwig_treated, metric='euclidean')
linkage_array = scipy.cluster.hierarchy.linkage(distance_matrix, method = 'ward')
optimised_order=scipy.cluster.hierarchy.optimal_leaf_ordering(linkage_array, distance_matrix, metric='euclidean')
new_order = scipy.cluster.hierarchy.leaves_list(optimised_order)
aligned_bigwig_overlay_sorted = []
for i in new_order:
    aligned_bigwig_overlay_sorted.append(aligned_bigwig_overlay[i])
aligned_bigwig_overlay_sorted = np.array(aligned_bigwig_overlay_sorted)
# %%
metadata_with_te_age_sorted=metadata_with_te_age.iloc[new_order]
# %%
bigwig_front_clustered = []
bigwig_back_clustered = []
for i in new_order:
    bigwig_front_clustered.append(bigwig_front_sorted[i])
    bigwig_back_clustered.append(bigwig_back_sorted[i])
#%%
fused_bigwig_mapped_sorted = []
for i in range(len(aligned_bigwig_overlay)):
    fused_bigwig_mapped_sorted.append(np.concatenate((bigwig_front_clustered[i], aligned_bigwig_overlay_sorted[i], bigwig_back_clustered[i])))
fused_bigwig_mapped_sorted = np.array(fused_bigwig_mapped_sorted)
#%%
age_colorcode = {0:'#bf6f93',6.7:'#79b743',9.06:'#b262cb',15.76:'#cea240',20.19:'#6566d1',29.44:'#738139',43.2:'#ca489e',73.8:'#51ad7b',76:'#cf4766',82:'#45b0cf',90:'#d25337',96:'#7981c7',105:'#b97446'}
subfam_colorcode={'MER11A': 'red', 'MER11B':'blue', 'MER11C':'yellow'}
metadata_with_te_age_sorted[['subfam', 'subfam_order']] = metadata_with_te_age_sorted['id'].str.split('_', expand=True)
row_color_subfam_sorted=metadata_with_te_age_sorted.subfam.map(subfam_colorcode)
row_color_age_sorted=metadata_with_te_age_sorted.te_age.map(age_colorcode)
row_colors_sorted=pd.DataFrame({'subfam':row_color_subfam_sorted,'te_age':row_color_age_sorted}).reset_index(drop=True)
#%%
graphical_object_sorted=sns.clustermap(pd.DataFrame(fused_bigwig_mapped_sorted),row_colors=row_colors_sorted, row_cluster=False, col_cluster=False, cmap =  "vlag_r",xticklabels =500, yticklabels = 500, annot = False, vmin=-1,vmax=1)
graphical_object_sorted.fig.subplots_adjust(left=0.05)
graphical_object_sorted.ax_cbar.set_position((0.1, .08, .02, .4))
graphical_object_sorted.ax_cbar.set_ylabel('phyloP')
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
graphical_object_sorted.ax_heatmap.set_title("MER11 phyloP overlay clustered")
plt.setp(graphical_object_sorted.ax_heatmap.set_xlabel("position (bp)"))
plt.setp(graphical_object_sorted.ax_heatmap.set_ylabel("sequences"))
plt.show()

# %%
metadata_sorted=metadata_with_te_age_sorted.reset_index(drop=True)
#%%
aligned_sorted = []
for i in new_order:
    aligned_sorted.append(aligned_filtered[i])
aligned_sorted = np.array(aligned_sorted)
# %%
metadata_sorted_mer11a=metadata_sorted[metadata_sorted.subfam == 'MER11A']
metadata_sorted_mer11b=metadata_sorted[metadata_sorted.subfam == 'MER11B']
metadata_sorted_mer11c=metadata_sorted[metadata_sorted.subfam == 'MER11C']
# %%
fused_bigwig_mer11a=fused_bigwig_mapped_sorted[metadata_sorted_mer11a.index]
aligned_sorted_mer11a=aligned_sorted[metadata_sorted_mer11a.index]
fused_bigwig_mer11b=fused_bigwig_mapped_sorted[metadata_sorted_mer11b.index]
aligned_sorted_mer11b=aligned_sorted[metadata_sorted_mer11b.index]
fused_bigwig_mer11c=fused_bigwig_mapped_sorted[metadata_sorted_mer11c.index]
aligned_sorted_mer11c=aligned_sorted[metadata_sorted_mer11c.index]
#%%
normalisation_mask_mer11a = np.count_nonzero(aligned_sorted_mer11a, axis=0)
border_count_mer11a = np.empty(500, dtype=np.int)
border_count_mer11a.fill(len(aligned_sorted_mer11a))
normaliser_mer11a=np.concatenate((border_count_mer11a,normalisation_mask_mer11a,border_count_mer11a))
normalisation_mask_mer11b = np.count_nonzero(aligned_sorted_mer11b, axis=0)
border_count_mer11b = np.empty(500, dtype=np.int)
border_count_mer11b.fill(len(aligned_sorted_mer11b))
normaliser_mer11b=np.concatenate((border_count_mer11b,normalisation_mask_mer11b,border_count_mer11b))
normalisation_mask_mer11c = np.count_nonzero(aligned_sorted_mer11c, axis=0)
border_count_mer11c = np.empty(500, dtype=np.int)
border_count_mer11c.fill(len(aligned_sorted_mer11c))
normaliser_mer11c=np.concatenate((border_count_mer11c,normalisation_mask_mer11c,border_count_mer11c))
#%%
fused_bigwig_treated_mer11a=np.nan_to_num(fused_bigwig_mer11a)
sum_fused_bigwig_mer11a = np.sum(fused_bigwig_treated_mer11a, axis = 0)
fused_bigwig_averaged_mer11a = sum_fused_bigwig_mer11a/normaliser_mer11a
fused_bigwig_treated_mer11b=np.nan_to_num(fused_bigwig_mer11b)
sum_fused_bigwig_mer11b = np.sum(fused_bigwig_treated_mer11b, axis = 0)
fused_bigwig_averaged_mer11b = sum_fused_bigwig_mer11b/normaliser_mer11b
fused_bigwig_treated_mer11c=np.nan_to_num(fused_bigwig_mer11c)
sum_fused_bigwig_mer11c = np.sum(fused_bigwig_treated_mer11c, axis = 0)
fused_bigwig_averaged_mer11c = sum_fused_bigwig_mer11c/normaliser_mer11c
#%%
fig, axs = plt.subplots(3,1,figsize=(10,6), sharex=True, sharey=True)
axs[0].fill_between(range(len(fused_bigwig_averaged_mer11a)), fused_bigwig_averaged_mer11a, color='grey')
axs[0].set_title('MER11A')
axs[1].fill_between(range(len(fused_bigwig_averaged_mer11b)), fused_bigwig_averaged_mer11b, color='grey')
axs[1].set_title('MER11B')
axs[2].fill_between(range(len(fused_bigwig_averaged_mer11c)), fused_bigwig_averaged_mer11c, color='grey')
axs[2].set_title('MER11C')
plt.setp(axs, ylim=(-1,1), xmargin=0)
#fig.margins(x=0, y=0)
#fig.set_ylim(0,1)
#fig.set_xlabel('position (bp)')
#fig.set_ylabel('normalised alt allele ratio')
#ig.set_title('MER11 bigwig overlay')
plt.show()

# %%
