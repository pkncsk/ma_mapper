# %%
from Bio import SeqIO
from Bio.AlignIO import MafIO
from collections import defaultdict
from Bio import AlignIO

# !/usr/bin/env python3
import os
import sys
sys.path.append('/home/pc575/rds/vcfcollector/variantcollector')
import importlib
import main
importlib.reload(main)
import numpy as np
import pandas as pd
### main
importlib.reload(main)
#%%
alignments = AlignIO.read('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/_housekeeping/soma_results/L1PA2/mafft_alignment.2024Mar12.fasta', "fasta")

#%%
alignment=main.alignment_init(
    output_dir = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/_housekeeping/soma_results/L1PA2/',
    source_bed = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/housekeeping/hg38_transposons_20140131_somatic.bed',
    target_group = [#'L1HS',
    'L1PA2',
    #'L1PA3',
    #'L1PA4',
    #'L1PA5',
    #'L1PA6',
    #'L1PA7',
    #'L1PA8',
    #'L1PA9',
    #'LTR10C',
    ],
    genome = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/hg38.fa',
    #sequence_number = 3000,
    group_level = 'subfamily', 
    force_overwrite = True,
    #randomise=True,
    #force_length=3000,
    #seed=65536,
    defragmentation= 50
)

#%%
##NOTE: MULTIPLE ALIGNMENT
alignment, metadata=main.parse_alignment(
    aligned_sequence_flie = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/_housekeeping/soma_results/L1PA2/mafft_alignment.2024Mar12.fasta',
    output_metadata_file = True
    )
#%%
alignment_sorted, metadata_sorted = main.sort_alignment_per_family(
    align_data = alignment,
    metadata = metadata,
    ascend= True,
)
filters = main.generate_filter(
    alignment_array = alignment_sorted
)
alignment_filtered=alignment_sorted[np.ix_(filters[0],filters[1])]
metadata_filtered=metadata_sorted.loc[filters[0]].copy(deep=True)
metadata_filtered['sort']=pd.Series(range(metadata_filtered.shape[0]), index=metadata_filtered.index)

main.alignment_plot(alignment_filtered)

# %%
#NOTE: PADDING 2000 bp +-
importlib.reload(main)
alignment_padded, metadata_padded, filter_padded=main.pad_data(
    sorted_alignment=alignment_sorted,
    sorted_metadata=metadata_sorted,
    filters=filters,
    padding_width=500
)
alignment_filtered=alignment_padded[np.ix_(filter_padded[0],filter_padded[1])]
metadata_filtered=metadata_padded.loc[filter_padded[0]].copy(deep=True)
metadata_filtered['sort']=pd.Series(range(metadata_filtered.shape[0]), index=metadata_filtered.index)

main.alignment_plot(alignment_filtered)
#%%
###NOTE: EVO
variant_data = main.process_variant_maf_signal(
    metadata_sorted = metadata_padded,
    coverage_count= False,
)
#%%
variant_data_filtered = main.apply_filtering(
    sorted_alignment_array = alignment_padded,
    signal_data = variant_data,
    filters = filter_padded
)
#%%
nonzero_row = np.count_nonzero(variant_data_filtered, axis=0)
#%%
processed_data = main.filtered_data_postprocessing(
    filtered_data = [variant_data_filtered],
    normalization_method = 'nonzero',
    ratio_calculation=False,
    #smoothing_method='gaussian',
    #smoothing_width=3 
)
#%%
main.ratio_plot(processed_data[0], 
    #minorticks_setting = 3, 
    smoothing= False, 
    color_list=['grey'],
    figsize = (24,6),
    ymax=1)
#%%
metadata.to_csv('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/l1pa2/old.txt', sep='\t', index= False)
#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import fetch_sequence
from ma_mapper import mafft_align
from ma_mapper import fetch_data
#%%
source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_fasta/hg38.fa'
metadata = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old.txt'
fasta_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/l1pa2_old.fasta'

fetch_sequence.fetch_sequence(metadata,source_fasta,output_filepath =fasta_file, save_to_file= True,custom_id= False)
#%%
fasta_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/_housekeeping/soma_results/L1PA2/sequence.2022Jun06.fasta'
mafft_align.mafft_wrapper(fasta_file, nthread = 40)
# %%
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import gcf
from ma_mapper import mapper
#input_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/_housekeeping/soma_results/L1PA2/sequence.2022Jun06.fasta.aligned'
input_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/_housekeeping/soma_results/L1PA2/mafft_alignment.2024Mar12.fasta'
aligned_parsed = mapper.parse_alignment(input_filepath, save_to_file= False)
#%%
records = mapper.load_alignment_file(input_filepath)
dump_list = []
chrom_list = []
start_list = []
end_list = []
strand_list = []
for record in records:
    subfamily, position_ = record.name.split(sep='::')
    chrom, coords=position_.split(sep=':')
    coord, strand_=coords.split(sep='(')
    start, end = coord.split(sep='-')
    strand = strand_.split(')')[0]
    dump_list.append(subfamily)
    chrom_list.append(chrom)
    start_list.append(int(start))
    end_list.append(int(end))
    strand_list.append(strand)
new_metadata_dict = {'chrom':chrom_list,'start':start_list,'end':end_list,'strand':strand_list,'id':dump_list}
metadata_aligned = pd.DataFrame(new_metadata_dict)
metadata_aligned['original_order'] = metadata_aligned.index
#%%
filters=mapper.create_filter(aligned_parsed)
row_filter = filters[0]
col_filter = filters[1]
aligned_filtered=aligned_parsed[np.ix_(row_filter,col_filter)]
metadata_aligned_filtered=metadata_aligned.iloc[row_filter,:]
#%%
map_color = ['grey','green','yellow','blue','red']
custom_cmap = LinearSegmentedColormap.from_list('Custom', map_color, len(map_color))
framed_alignment=pd.DataFrame(aligned_filtered)
graphical_object=sns.clustermap(framed_alignment, row_cluster=False, col_cluster=False, cmap =  custom_cmap,xticklabels =aligned_filtered.shape[1]-1, yticklabels = 500, annot = False)
graphical_object.ax_row_dendrogram.set_visible(False)
graphical_object.ax_col_dendrogram.set_visible(False)
graphical_object.fig.subplots_adjust(left=0.05)
graphical_object.ax_cbar.set_position((0.15, .08, .02, .4))
graphical_object.ax_cbar.set_ylabel('allele')
graphical_object.ax_cbar.set_yticklabels(['gap','','A','','C','','T','','G'])
graphical_object.ax_heatmap.set_title("l1pa2 alignment")

plt.setp(graphical_object.ax_heatmap.set_xlabel("position (bp)"))
plt.setp(graphical_object.ax_heatmap.set_ylabel("sequences"))
plt.gcf().set_size_inches(24, 6)
plt.show()
#%% find border for each id
# make new metadata 
temp_dict = {'chrom':metadata_aligned_filtered.chrom.to_list(),
             'start':metadata_aligned_filtered.start.to_list(),'end':metadata_aligned_filtered.start.to_list(),
             'strand':metadata_aligned_filtered.strand.to_list(),
             'id':metadata_aligned_filtered.original_order.to_list()}
low_border_metadata = pd.DataFrame(temp_dict)
low_border_metadata.start = low_border_metadata.start-500
# %%
source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/hg38.fa'
low_border_records=fetch_sequence.fetch_sequence(low_border_metadata,source_fasta, custom_id= False)

# %%
temp_dict = {'chrom':metadata_aligned_filtered.chrom.to_list(),
             'start':metadata_aligned_filtered.end.to_list(),
             'end':metadata_aligned_filtered.end.to_list(),
             'strand':metadata_aligned_filtered.strand.to_list(),
             'id':metadata_aligned_filtered.original_order.to_list()}
high_border_metadata = pd.DataFrame(temp_dict)
high_border_metadata.end = high_border_metadata.end+500
#%%
source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/_housekeeping/data/hg38.fa'
high_border_records=fetch_sequence.fetch_sequence(high_border_metadata,source_fasta, custom_id= False)
# %%
front_list = []
back_list = []
for idx, strand in enumerate(metadata_aligned_filtered.strand.to_list()):
    if strand == '+':
        front_list.append(low_border_records[idx])
        back_list.append(high_border_records[idx])
    else:
        front_list.append(high_border_records[idx])
        back_list.append(low_border_records[idx])
    
# %%
front_parsed = mapper.parse_alignment(front_list, save_to_file= False)
back_parsed = mapper.parse_alignment(back_list, save_to_file= False)
# %%
fused_parsed = []
for i in range(len(aligned_filtered)):
    fused_parsed.append(np.concatenate((front_parsed[i], aligned_filtered[i], back_parsed[i])))
fused_parsed = np.array(fused_parsed)
#%%
map_color = ['grey','green','yellow','blue','red']
custom_cmap = LinearSegmentedColormap.from_list('Custom', map_color, len(map_color))
framed_alignment=pd.DataFrame(fused_parsed)
graphical_object=sns.clustermap(framed_alignment, row_cluster=False, col_cluster=False, cmap =  custom_cmap,xticklabels =aligned_filtered.shape[1]-1, yticklabels = 500, annot = False)
graphical_object.ax_row_dendrogram.set_visible(False)
graphical_object.ax_col_dendrogram.set_visible(False)
graphical_object.fig.subplots_adjust(left=0.05)
graphical_object.ax_cbar.set_position((0.15, .08, .02, .4))
graphical_object.ax_cbar.set_ylabel('allele')
graphical_object.ax_cbar.set_yticklabels(['gap','','A','','C','','T','','G'])
graphical_object.ax_heatmap.set_title("l1pa2 alignment")

plt.setp(graphical_object.ax_heatmap.set_xlabel("position (bp)"))
plt.setp(graphical_object.ax_heatmap.set_ylabel("sequences"))
plt.gcf().set_size_inches(24, 6)
plt.show()
#%% extract maf
maf_mapped=fetch_data.fetch_maf(metadata_input= metadata_aligned, maf_input='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/241genomes/241-mammalian-2020v2b.maf', separated_maf = True,target_species = 'Homo_sapiens', custom_id= False)
#%%
aligned_maf_overlay=mapper.map_data(maf_mapped, aligned_parsed, filters = filters)
metadata_aligned_filtered=metadata_aligned.iloc[row_filter,:]
# %%
low_border_maf_mapped=fetch_data.fetch_maf(metadata_input= low_border_metadata, maf_input='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/241genomes/241-mammalian-2020v2b.maf', separated_maf = True,target_species = 'Homo_sapiens', custom_id= True)
# %%
high_border_maf_mapped=fetch_data.fetch_maf(metadata_input= high_border_metadata, maf_input='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/241genomes/241-mammalian-2020v2b.maf', separated_maf = True,target_species = 'Homo_sapiens', custom_id= True)
#%%
maf_front_list = []
maf_back_list = []
for idx, strand in enumerate(metadata_aligned_filtered.strand.to_list()):
    if strand == '+':
        maf_front_list.append(low_border_maf_mapped[idx])
        maf_back_list.append(high_border_maf_mapped[idx])
    else:
        maf_front_list.append(high_border_maf_mapped[idx])
        maf_back_list.append(low_border_maf_mapped[idx])
#%%
fused_maf_mapped = []
for i in range(len(aligned_maf_overlay)):
    fused_maf_mapped.append(np.concatenate((maf_front_list[i], aligned_maf_overlay[i], maf_back_list[i])))
fused_maf_mapped = np.array(fused_maf_mapped)
#%%
graphical_object_sorted=sns.clustermap(pd.DataFrame(fused_maf_mapped), row_cluster=False, col_cluster=False, cmap =  "viridis",xticklabels =500, yticklabels = 500, annot = False)
graphical_object_sorted.fig.subplots_adjust(left=0.05)
graphical_object_sorted.ax_cbar.set_position((0.15, .08, .02, .4))
graphical_object_sorted.ax_cbar.set_ylabel('normalised_alt_allele')
col = graphical_object_sorted.ax_col_dendrogram.get_position()
graphical_object_sorted.ax_col_dendrogram.set_position([col.x0, col.y0, col.width*0.25, col.height*0.25])
#graphical_object.ax_row_dendrogram.set_visible(False)
#graphical_object.ax_col_dendrogram.set_visible(False)

graphical_object_sorted.ax_heatmap.set_title("l1pa2 MAF overlay")
plt.setp(graphical_object_sorted.ax_heatmap.set_xlabel("position (bp)"))
plt.setp(graphical_object_sorted.ax_heatmap.set_ylabel("sequences"))
plt.show()
#%%
aligned_filtered=aligned_parsed[np.ix_(row_filter,col_filter)]
normalisation_mask = np.count_nonzero(aligned_filtered, axis=0)
border_count = np.empty(500, dtype=np.int)
border_count.fill(len(maf_front_list))
normaliser=np.concatenate((border_count,normalisation_mask,border_count))
fused_maf_treated=np.nan_to_num(fused_maf_mapped)
sum_fused_maf_mapped = np.sum(fused_maf_treated, axis = 0)
fused_maf_averaged = sum_fused_maf_mapped/normaliser
# %%
plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600
fig, ax = plt.subplots(figsize=(10,3))
ax.fill_between(range(len(fused_maf_averaged)), fused_maf_averaged, color='grey')
ax.margins(x=0, y=0)
ax.set_ylim(0,1)
ax.set_xlabel('position (bp)')
ax.set_ylabel('averaged alt allele ratio')
ax.set_title('l1pa2 MAF overlay')
plt.show()
#%%
sum_maf = np.sum(fused_maf_mapped, axis =0)
normalizer=np.count_nonzero(fused_maf_mapped, axis =0)
fused_maf_normalized= sum_maf/normalizer
#%%
np.sum(variant_data_filtered[:,500:-500], axis=0)
np.sum(fused_maf_mapped[:,500:-500], axis=0)
# %%
plt.rcParams['figure.dpi'] = 75
plt.rcParams['savefig.dpi'] = 75
fig, ax = plt.subplots(figsize=(24,6))
ax.fill_between(range(len(fused_maf_normalized)), fused_maf_normalized, color='grey')
ax.margins(x=0, y=0)
ax.set_ylim(0,1)
ax.set_xlabel('position (bp)')
ax.set_ylabel('averaged alt allele ratio')
ax.set_title('l1pa2 MAF overlay')
plt.show()
# %%
import pybedtools
metadata_for_bed=metadata.iloc[:,1:4]
metadata_for_bed['name'] = metadata.family
metadata_for_bed['score'] = 0
metadata_for_bed['strand'] = metadata.strand
metadata_bed=pybedtools.BedTool.from_dataframe(metadata_for_bed)
# %%
new_metadata_df=pd.read_csv('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/l1pa2_coord_with_id.txt',sep='\t')
# %%
metadata_for_bed=new_metadata_df.iloc[:,0:3]
metadata_for_bed['name'] = new_metadata_df.id
metadata_for_bed['score'] = 0
metadata_for_bed['strand'] = new_metadata_df.strand
# %%
new_metadata_bed=pybedtools.BedTool.from_dataframe(metadata_for_bed)
# %%
temp_intersect=metadata_bed.intersect(new_metadata_bed,loj=True, wa = True, wb =True, s=True)
intersect_df=temp_intersect.to_dataframe()
# %%
main.alignment_plot(aligned_filtered)
# %%
#%%
###NOTE: EVO
variant_data_test = main.process_variant_maf_signal(
    metadata_sorted = metadata_sorted,
    coverage_count= False,
)
# %%
variant_data_test
# %%
maf_mapped
#%%

#def extract_variant_freq(chrom, start, end, family, sort, strand):
target_species = 'Homo_sapiens'
chrom = 'chr14'
start = 37350446
end = 37350456
strand = '-'
coverage_count = False
filepath = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/241genomes/241-mammalian-2020v2b.maf.'+str(chrom)
id_acc = target_species+'.'+chrom
idx = MafIO.MafIndex(filepath+".mafindex", filepath, id_acc)
if strand=='-':
    strand = -1
else:
    strand = 1
#offseting the surro
# unding sequences
#start = bed_record.start-2500
#end = bed_record.end+2500
#results =idx.get_spliced([start],[end],strand)
results =idx.get_spliced([start],[end],strand)
#NOTE: the get_spliced already take care of the minus strand 
#NOTE: data can be used without additional reversing
test_list = []
collector = defaultdict()
for seqrec in results:
    #print(seqrec.id.split('.')[0])
    if seqrec.id.split('.')[0] == target_species:
        #print(seqrec.seq)
        ref_seq = seqrec.seq
        #test_list.append(np.array(list(ref_seq)))
        collector[seqrec.id.split('.')[0]]= np.array(list(seqrec.seq))
        for seqrec in results:
            if seqrec.id.split('.')[0] != target_species:
                #test_list.append(np.array(list(seqrec.seq)))
                print(seqrec.id,seqrec.seq)
                collector[seqrec.id.split('.')[0]]= np.array(list(seqrec.seq))
#testarray_transposed=np.array(test_list).transpose()
testarray_transposed=np.array(list(collector.values())).transpose()
alt_freq_array=[]
for ref_pos in testarray_transposed:
    (unique, counts) = np.unique(ref_pos, return_counts=True)
    frequencies = dict(zip(unique, counts))
    ref_allele=ref_pos[0]
    total = 0
    alt_count  = 0
    for key in frequencies:
        if key != '-':
            total = total+frequencies[key]
            if key != ref_allele:
                alt_count =alt_count+frequencies[key]
    if coverage_count is True:
        alt_freq = total
    else:
        alt_freq=alt_count/total
    #print(ref_allele+':'+str(frequencies[ref_allele])+' alt:'+str(alt_count)+' total:'+str(total))
    alt_freq_array.append(alt_freq)
print(alt_freq_array)

# %%
target_species = 'Homo_sapiens'
chrom = 'chr14'
start = 37350446
end = 37350456
strand = '-'
coverage_count = False
age_arg = None
age_table_file = None
maf_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/241genomes/241-mammalian-2020v2b.maf.chr14'
if age_arg is not None:
    if isinstance(age_table_file, str):
        if (os.path.isfile(age_table_file) == True):
            species = pd.read_table(age_table_file, sep='\t')
        else:
            print('metadata file not found')
    else:
        species = age_table_file
maf_id = target_species+'.'+chrom
from Bio.AlignIO import MafIO
index_maf = MafIO.MafIndex(maf_file+".mafindex", maf_file, maf_id) 
if strand =='-':
    n_strand = -1
    results =index_maf.get_spliced([start],[end],n_strand)
else:
    n_strand = 1
    results =index_maf.get_spliced([start],[end],n_strand)
collector = {}
for seqrec in results:
    if seqrec.id.split('.')[0] == target_species:  
        collector[seqrec.id.split('.')[0]]= np.array(list(seqrec.seq.lower()))
        age_measure = 0
        for seqrec in results:
            if seqrec.id.split('.')[0] != target_species:
                if age_arg is None:
                    print(seqrec.id,seqrec.seq)
                    collector[seqrec.id.split('.')[0]]= np.array(list(seqrec.seq.lower()))
                elif age_arg == 'calibrate':
                    species_age = species[species.meta_name == seqrec.id.split('.')[0]]['Estimated Time (MYA)'].to_list()[0]
                    if age >= species_age:
                        collector[seqrec.id.split('.')[0]]= np.array(list(seqrec.seq.lower()))
                elif age_arg == 'extract':
                    species_age = species[species.meta_name == seqrec.id.split('.')[0]]['Estimated Time (MYA)'].to_list()[0]
                    if species_age >= age_measure:
                        age_measure = species_age

array_transposed=np.array(list(collector.values())).transpose()
alt_freq_array=[]
for ref_pos in array_transposed:
    (unique, counts) = np.unique(ref_pos, return_counts=True)
    frequencies = dict(zip(unique, counts))
    ref_allele=ref_pos[0]
    total = 0
    alt_count  = 0
    for key in frequencies:
        if key != '-':
            total = total+frequencies[key]
            if key != ref_allele:
                alt_count =alt_count+frequencies[key]
    if coverage_count is True:
        alt_freq = total
    else:
        alt_freq=alt_count/total
    alt_freq_array.append(alt_freq)
print(alt_freq_array)
# %%
