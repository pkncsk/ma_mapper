#%%
import sys
import pandas as pd
import numpy as np
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import mapper
from ma_mapper import fetch_data
import pybedtools
#%%
input_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11.fasta.aligned'
aligned_parsed = mapper.parse_alignment(input_filepath, save_to_file= False)
metadata_aligned = mapper.extract_metadata_from_alignment(input_filepath)
#%%
metadata_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11_coord_with_id.txt'
#%%
metadata_df = pd.read_csv(metadata_filepath, sep='\t')
original_order = metadata_df.iloc[:,4].unique()
#%%
meta_id=metadata_df.id.unique()
metadata_for_bed=metadata_df.iloc[:,0:3]
metadata_for_bed['name'] = metadata_df.id
metadata_for_bed['score'] = 10
metadata_for_bed['strand'] = metadata_df.strand
# %%
bed_filepath ='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/tfbs_homer/Gata6(Zf).bed'
metadata_bed=pybedtools.BedTool.from_dataframe(metadata_for_bed)
bed_file = pybedtools.BedTool(bed_filepath)
# %%
temp_intersect=metadata_bed.intersect(bed_file,loj=True, wa = True, wb =True, s=True)
intersect_df=temp_intersect.to_dataframe()
#%%
bed_out = []
for unique_id in meta_id:
    intersect_by_id=intersect_df[intersect_df.name == unique_id]
    unique_start=intersect_by_id.start.unique()
    strand=intersect_by_id.strand.unique()[0]
    drawings = []
    for start in unique_start:
        intersect_by_id_by_start=intersect_by_id[intersect_by_id.start == start]
        canvas_start = intersect_by_id_by_start.start.unique()[0]
        canvas_end = intersect_by_id_by_start.end.unique()[0]
        if(len(intersect_by_id_by_start.end.unique())>1):
            print('multiple length frags')
        canvas = np.zeros(canvas_end-canvas_start)
        for idx,row in intersect_by_id_by_start.iterrows():
            if int(row.thickEnd) != -1:
                intersect_start = int(row.thickEnd)
                intersect_end = int(row.itemRgb)
                intersect_score = int(row.blockSizes)
                #capping out of bound position
                if intersect_end > canvas_end:
                    intersect_end = canvas_end
                if intersect_start < canvas_start:
                    intersect_start = canvas_start
                if strand == '+':
                    start_intersect_position = intersect_start - canvas_start
                    end_intersect_postion = intersect_end - canvas_start

                else:
                    start_intersect_position = canvas_end - intersect_end
                    end_intersect_postion = canvas_end - intersect_start
                canvas[start_intersect_position:end_intersect_postion] += intersect_score
        #if strand == '-':
            #canvas = np.flip(canvas)
        drawings.append(canvas)
    if strand == '+':
        bed_out.append(np.concatenate(drawings))
    else:
        reverse_order=np.flip(drawings, 0)
        bed_out.append(np.concatenate(reverse_order))
#%%
bed_mapped_sorted = []
for idx, row in metadata_aligned.iterrows():
    bed_mapped_sorted.append(bed_out[np.where(original_order == row.id)[0][0]])
#%%
filters=mapper.create_filter(aligned_parsed,row_threshold=0.1, col_threshold=0.1,col_content_threshold=0.1)
#%%
aligned_bed_overlay=mapper.map_data(bed_mapped_sorted, aligned_parsed, filters = filters)
# %%
row_filter = filters[0]
col_filter = filters[1]
aligned_filtered=aligned_parsed[np.ix_(row_filter,col_filter)]
metadata_aligned_filtered=metadata_aligned.iloc[row_filter,:]

# %%
overlap_id=intersect_df[intersect_df.blockSizes != -1].name.unique()
# %%
overlap_index=metadata_aligned_filtered[metadata_aligned_filtered.id.isin(overlap_id)].reset_index().index
nonoverlap_index=metadata_aligned_filtered[~metadata_aligned_filtered.id.isin(overlap_id)].reset_index().index
#%%
import sys
import pandas as pd
import numpy as np
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import mapper
from ma_mapper import fetch_data
import matplotlib.pyplot as plt
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
bigwig_mapped=fetch_data.fetch_bigwig(metadata_input= metadata_filepath, bigwig_input='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_bigwig/241-mammalian-2020v2.bigWig', custom_id= True)
#%%
metadata_df = pd.read_csv(metadata_filepath, sep='\t')
original_order = metadata_df.iloc[:,4].unique()
#%%
bigwig_mapped_sorted = []
for idx, row in metadata_aligned.iterrows():
    bigwig_mapped_sorted.append(bigwig_mapped[np.where(original_order == row.id)[0][0]])
#%%
filters=mapper.create_filter(aligned_parsed,row_threshold=0.1, col_threshold=0.1,col_content_threshold=0.1)
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
overlap_bigwig=fused_bigwig_mapped[overlap_index,:]
nonoverlap_bigwig=fused_bigwig_mapped[nonoverlap_index,:]
# %%
from scipy import stats
#%%
stat_v, p_value=stats.ttest_ind(overlap_bigwig,nonoverlap_bigwig, axis =0,nan_policy='omit')
# %%S
#plt.rcParams['figure.dpi'] = 600
#plt.rcParams['savefig.dpi'] = 600
#fig, ax = plt.subplots(figsize=(10,3))
#ax.fill_between(range(len(p_value)), -np.log10(p_value), color='grey')
#ax.margins(x=0, y=0)
#ax.set_ylim(-1,1)S
#ax.set_xlabel('position (bp)')
#ax.set_ylabel('normalised phyloP')
#ax.set_title('MER11 phyloP overlay')
#plt.show()
# %%
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import gcf
#%%
#make transparent cmap
ncolors = 256
color_array = plt.get_cmap('Blues')(range(ncolors))
# change alpha values
color_array[:,-1] = np.linspace(0,1.0,ncolors)
# create a colormap object
blues_alpha = LinearSegmentedColormap.from_list(name='rainbow_alpha',colors=color_array)
# %%
#fig = plt.figure(figsize=(10,12))
#grid = fig.add_gridspec(nrows = 6, ncols = 1, hspace=0)
#ax0 = fig.add_subplot(grid[0:5,:])
#ax1 = fig.add_subplot(grid[5,:])
#square_aspect = aligned_bigwig_overlay.shape[1]/aligned_bigwig_overlay.shape[0]
#heatmap_core = ax0.imshow(aligned_bigwig_overlay, aspect = 'auto', cmap='RdBu', vmin = -1, vmax=1)
#ax0.imshow(aligned_bed_overlay, aspect = 'auto', cmap = blues_alpha, vmax = 5)
#ax1.fill_between(range(len(p_value[500:-500])), -np.log10(p_value[500:-500]), color='grey')
#ax1.margins(x=0, y=0)
#plt.show()
# %%
overlap_bigwig[np.random.choice(overlap_bigwig.shape[0], nonoverlap_bigwig.shape[0], replace=False), :]
# %%
time_permute = 100
p_value_permuted = []
for i in range(time_permute):
    print(i)
    subsampled_overlap_bigwig=overlap_bigwig[np.random.choice(overlap_bigwig.shape[0], nonoverlap_bigwig.shape[0], replace=False), :]
    stat_v, p_value=stats.ttest_ind(subsampled_overlap_bigwig,nonoverlap_bigwig, axis =0,nan_policy='omit')
    p_value_permuted.append(p_value)
# %%
p_val_np=np.array(p_value_permuted)
p_val_permulted_mean=np.mean(p_value_permuted, axis = 0)
# %%
# %%
#fig = plt.figure(figsize=(10,12))
#grid = fig.add_gridspec(nrows = 6, ncols = 1, hspace=0)
#ax0 = fig.add_subplot(grid[0:5,:])
#ax1 = fig.add_subplot(grid[5,:])
#square_aspect = aligned_bigwig_overlay.shape[1]/aligned_bigwig_overlay.shape[0]
#heatmap_core = ax0.imshow(aligned_bigwig_overlay, aspect = 'auto', cmap='RdBu', vmin = -1, vmax=1)
#ax0.imshow(aligned_bed_overlay, aspect = 'auto', cmap = blues_alpha, vmax = 5)
#ax1.fill_between(range(len(p_val_permulted_mean[500:-500])), -np.log10(p_val_permulted_mean[500:-500]), color='grey')
#ax1.margins(x=0, y=0)
#plt.show()

# %%
#overlap_bigwig_trimmed=overlap_bigwig[:,500:-500]
#overlapped_bed=aligned_bed_overlay[overlap_index,:]
#fig = plt.figure(figsize=(10,12))
#grid = fig.add_gridspec(nrows = 6, ncols = 1, hspace=0)
#ax0 = fig.add_subplot(grid[0:5,:])
#ax1 = fig.add_subplot(grid[5,:])
#square_aspect = aligned_bigwig_overlay.shape[1]/aligned_bigwig_overlay.shape[0]
#heatmap_core = ax0.imshow(overlap_bigwig_trimmed, aspect = 'auto', cmap='RdBu', vmin = -1, vmax=1)
#ax0.imshow(overlapped_bed, aspect = 'auto', cmap = blues_alpha, vmax = 5)
#ax1.fill_between(range(overlap_bigwig_trimmed.shape[1]), np.mean(overlap_bigwig_trimmed, axis = 0), color='grey')
#ax1.margins(x=0, y=0)
#plt.show()
#%%
np.sum(aligned_bed_overlay, axis = 0)
#%%
from scipy.signal import find_peaks
peaks, peaks_info = find_peaks(np.sum(aligned_bed_overlay, axis = 0), height=200, width =7)

plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600
fig, ax = plt.subplots(figsize=(10,3))
ax.plot(peaks, np.sum(aligned_bed_overlay, axis = 0)[peaks], 'x')
ax.fill_between(range(aligned_bed_overlay.shape[1]), np.sum(aligned_bed_overlay, axis = 0), color='grey')
ax.margins(x=0, y=0)
#ax.set_ylim(-1,1)S
ax.set_xlabel('position (bp)')
ax.set_ylabel('normalised phyloP')
ax.set_title('MER11 phyloP overlay')
plt.show()
#%% find peak overlap
peaks = np.array([248, 790, 827])
#peak_idx = 755
p_values = []
for peak_idx in peaks:
    print(peak_idx)
    overlap_idx = []
    nonoverlap_idx = []
    for idx in range(len(aligned_bed_overlay)):
        if aligned_bed_overlay[idx][peak_idx]> 0:
            overlap_idx.append(idx)
        else:
            nonoverlap_idx.append(idx)
    overlap_bigwig=fused_bigwig_mapped[overlap_idx,:]
    nonoverlap_bigwig=fused_bigwig_mapped[nonoverlap_idx,:]
    stat_v, p_value=stats.ttest_ind(overlap_bigwig,nonoverlap_bigwig, axis =0,nan_policy='omit')
    p_values.append(p_value)
#%%
#overlap_bigwig_trimmed=overlap_bigwig[:,500:-500]
pad = np.zeros((aligned_bigwig_overlay.shape[0],500))
#overlapped_bed_overlapped=aligned_bed_overlay[overlap_idx,:]
aligned_bed_padded = np.concatenate((pad, aligned_bed_overlay, pad), axis = 1)
#%%
plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600
fig = plt.figure(figsize=(10,14))
grid = fig.add_gridspec(nrows = 9, ncols = 1, hspace=0)
ax0 = fig.add_subplot(grid[0:5,:])
ax1 = fig.add_subplot(grid[5,:])
ax2 = fig.add_subplot(grid[6,:])
ax3 = fig.add_subplot(grid[7,:])
ax4 = fig.add_subplot(grid[8,:])
#square_aspect = aligned_bigwig_overlay.shape[1]/aligned_bigwig_overlay.shape[0]
heatmap_core = ax0.imshow(fused_bigwig_mapped, aspect = 'auto', cmap='RdBu', vmin = -1, vmax=1)
ax0.imshow(aligned_bed_padded, aspect = 'auto', cmap = blues_alpha, vmax = 5)
ax1.fill_between(range(aligned_bed_padded.shape[1]), np.sum(aligned_bed_padded, axis = 0), color='grey')
ax1.margins(x=0, y=0)
ax2.fill_between(range(aligned_bed_padded.shape[1]), -np.log10(p_values[0]), color='grey')
ax2.axvline(x=500+248)
ax2.margins(x=0, y=0)
ax3.fill_between(range(aligned_bed_padded.shape[1]), -np.log10(p_values[1]), color='grey')
ax3.axvline(x=500+790)
ax3.margins(x=0, y=0)
ax4.fill_between(range(aligned_bed_padded.shape[1]), -np.log10(p_values[2]), color='grey')
ax4.axvline(x=500+827)
ax4.margins(x=0, y=0)
plt.show()
# %%
peak_idx = 721
overlap_idx = []
nonoverlap_idx = []
for idx in range(len(aligned_bed_overlay)):
    if aligned_bed_overlay[idx][peak_idx]> 0:
        overlap_idx.append(idx)
    else:
        nonoverlap_idx.append(idx)
overlap_bigwig=fused_bigwig_mapped[overlap_idx,:]
nonoverlap_bigwig=fused_bigwig_mapped[nonoverlap_idx,:]
#=np.nan_to_num(overlap_bigwig)
#nonoverlap_bigwig_treated=np.nan_to_num(nonoverlap_bigwig)
stat_v, p_value=stats.ttest_ind(overlap_bigwig,nonoverlap_bigwig, axis =0,nan_policy='omit')
#%%
plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600
fig, ax = plt.subplots(figsize=(10,3))
ax.fill_between(range(len(p_value)), -np.log10(p_value), color='grey')
ax.margins(x=0, y=0)
ax.set_xlabel('position (bp)')
ax.set_ylabel('p-value (-log10)')
ax.set_title('MER11 phyloP t-test')
#plt.show()
# %%
for i in range(len(p_value)):
    print(p_value[i]-p_values[1][i])
# %%
test = {'a':5,'b':4,'c':3}
# %%
