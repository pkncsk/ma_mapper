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
filters=mapper.create_filter(aligned_parsed)
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
#peaks = np.array([244, 721, 755])
peaks = [755]
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

# %%
#reversh hashing
hash_table = {0:'-',1:'a',2:'c',3:'t',4:'g'}
reharsh = np.vectorize(lambda x: hash_table[x])
aligned_filtered_str = reharsh(aligned_filtered)
alinged_transposed=aligned_filtered_str.transpose()
aligned_concensus = []
for ref_pos in alinged_transposed:
    (unique, counts) = np.unique(ref_pos, return_counts=True)
    frequencies = dict(zip(unique, counts))
    #print(frequencies)
    frequencies.pop('-', None)
    common_allele = max(frequencies, key=frequencies.get)
    aligned_concensus.append(common_allele)

#%%
overlap_aligned=aligned_filtered_str[overlap_idx,:]
nonoverlap_aligned=aligned_filtered_str[nonoverlap_idx,:]
overlap_transposed=overlap_aligned.transpose()
nonoverlap_transposed=nonoverlap_aligned.transpose()
overlap_count = []
nonoverlap_count = []
for i in range(len(aligned_concensus)): 
    ref_pos=overlap_transposed[i]
    (unique, counts) = np.unique(ref_pos, return_counts=True)
    frequencies = dict(zip(unique, counts))
    frequencies.pop('-', None)
    #print(i, frequencies)
    #total = sum(frequencies.values())
    #concensus_count = frequencies.get(aligned_concensus[i],0)
    #alt_count = total - concensus_count
    #overlap_count.append([concensus_count, alt_count])
    allele_count = []
    for key in ['a','c','t','g']:
        base_count = frequencies.get(key,0)
        allele_count.append(base_count)
    overlap_count.append(allele_count)
for i in range(len(aligned_concensus)):
    ref_pos=nonoverlap_transposed[i]
    (unique, counts) = np.unique(ref_pos, return_counts=True)
    frequencies = dict(zip(unique, counts))
    frequencies.pop('-', None)
    #total = sum(frequencies.values())
    #concensus_count = frequencies.get(aligned_concensus[i],0)
    #alt_count = total - concensus_count
    #nonoverlap_count.append([concensus_count, alt_count])
    allele_count = []
    for key in ['a','c','t','g']:
        base_count = frequencies.get(key,0)
        allele_count.append(base_count)
    nonoverlap_count.append(allele_count)    
# %%
from scipy.stats import chi2_contingency
p_val_chi2 = []
for i in range(len(aligned_concensus)):
    contingency_table = np.array([overlap_count[i],nonoverlap_count[i]])
    zero_columns = np.all(contingency_table == 0, axis=0)
    new_contingency_table = contingency_table[:, ~zero_columns]
    chi, p, dof, expected = chi2_contingency(new_contingency_table)
    p_val_chi2.append(p)
    #else: p_val_chi2.append(1)
# %%
plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600
fig, ax = plt.subplots(figsize=(10,3))
ax.fill_between(range(len(p_val_chi2)), -np.log10(p_val_chi2), color='grey')
ax.margins(x=0, y=0)
#ax.set_ylim(-1,1)S
ax.set_xlabel('position (bp)')
ax.set_ylabel('normalised phyloP')
ax.set_title('MER11 phyloP overlay')
plt.show()
# %%
