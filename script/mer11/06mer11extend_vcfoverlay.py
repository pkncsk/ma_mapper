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
vcf_mapped=fetch_data.fetch_vcf(metadata_input= metadata_filepath, vcf_input='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/vcf-gnomad/', vcf_format = 'gnomad', custom_id= True)
#%%
metadata_df = pd.read_csv(metadata_filepath, sep='\t')
original_order = metadata_df.iloc[:,4].unique()
#%%
vcf_mapped_sorted = []
for idx, row in metadata_aligned.iterrows():
    vcf_mapped_sorted.append(vcf_mapped[np.where(original_order == row.id)[0][0]])
#%%
filters=mapper.create_filter(aligned_parsed)
row_filter = filters[0]
col_filter = filters[1]
# %%
aligned_vcf_overlay=mapper.map_data(vcf_mapped_sorted, aligned_parsed, filters = filters)
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
low_border_vcf_mapped=fetch_data.fetch_vcf(metadata_input= low_border_metadata, vcf_input='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/vcf-gnomad/', vcf_format = 'gnomad', custom_id= True)
# %%
high_border_vcf_mapped=fetch_data.fetch_vcf(metadata_input= high_border_metadata, vcf_input='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/vcf-gnomad/', vcf_format = 'gnomad', custom_id= True)
#%%
vcf_front_list = []
vcf_back_list = []
for idx, strand in enumerate(strand_list):
    if strand == '+':
        vcf_front_list.append(low_border_vcf_mapped[idx])
        vcf_back_list.append(high_border_vcf_mapped[idx])
    else:
        vcf_front_list.append(high_border_vcf_mapped[idx])
        vcf_back_list.append(low_border_vcf_mapped[idx])
vcf_front_list = np.array(vcf_front_list, dtype = 'float')
vcf_back_list=np.array(vcf_back_list, dtype = 'float')
#%%
vcf_front_sorted = []
vcf_back_sorted = []
for idx, row in metadata_aligned.iterrows():
    if row_filter[idx]:
        vcf_front_sorted.append(vcf_front_list[np.where(original_order == row.id)[0][0]])
        vcf_back_sorted.append(vcf_back_list[np.where(original_order == row.id)[0][0]])
#%%
fused_vcf_mapped = []
for i in range(len(aligned_vcf_overlay)):
    fused_vcf_mapped.append(np.concatenate((vcf_front_sorted[i], aligned_vcf_overlay[i], vcf_back_sorted[i])))
fused_vcf_mapped = np.array(fused_vcf_mapped, dtype = 'float')
# %%
# %%
age_colorcode = {0:'#bf6f93',6.7:'#79b743',9.06:'#b262cb',15.76:'#cea240',20.19:'#6566d1',29.44:'#738139',43.2:'#ca489e',73.8:'#51ad7b',76:'#cf4766',82:'#45b0cf',90:'#d25337',96:'#7981c7',105:'#b97446'}
subfam_colorcode={'MER11A': 'red', 'MER11B':'blue', 'MER11C':'yellow'}
metadata_with_te_age[['subfam', 'subfam_order']] = metadata_with_te_age['id'].str.split('_', expand=True)
row_color_subfam_sorted=metadata_with_te_age.subfam.map(subfam_colorcode)
row_color_age_sorted=metadata_with_te_age.te_age.map(age_colorcode)
row_colors_sorted=pd.DataFrame({'subfam':row_color_subfam_sorted,'te_age':row_color_age_sorted}).reset_index(drop=True)

#%%
graphical_object_sorted=sns.clustermap(pd.DataFrame(fused_vcf_mapped),row_colors=row_colors_sorted, 
row_cluster=False, col_cluster=False, cmap =  "Blues",xticklabels =500, yticklabels = 500, annot = False, vmax=0.0005)
graphical_object_sorted.fig.subplots_adjust(left=0.05)
graphical_object_sorted.ax_cbar.set_position((0.1, .08, .02, .4))
graphical_object_sorted.ax_cbar.set_ylabel('normalised_alt_allele_freq')
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
graphical_object_sorted.ax_heatmap.set_title("MER11 vcf overlay")
plt.setp(graphical_object_sorted.ax_heatmap.set_xlabel("position (bp)"))
plt.setp(graphical_object_sorted.ax_heatmap.set_ylabel("sequences"))
plt.show()
#%%