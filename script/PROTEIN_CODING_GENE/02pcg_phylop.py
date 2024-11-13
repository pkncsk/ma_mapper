#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import extract_bigwig
from ma_mapper import mapper
#%%
subfamily = ['pcg_sliced']
alignment_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/alignment/'+subfamily[0]+'.fasta'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered)
# %%
coord_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/coord_internal_id/'+subfamily[0]+'.txt'
bigwig_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_bigwig/241-mammalian-2020v2.bigWig'

#phylop=extract_bigwig.fetch_bigwig(coord_file, bigwig_file, custom_id = True)
#%%
phylop=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig', extension_length=500)

# %%
from ma_mapper._unused import plot
plot.overlay_plot([alignment_filtered,phylop], metadata= metadata_age,data_vmin=[-1.0],data_vmax=[1.0], nucleotide_color='white',data_cmap=['RdBu'], plot_title='pcg_sliced:Zoonomia_phyloP', show_data_legend=True,data_label=['phyloP'], age_annotation=False, xlim=[400,550])

#%%
import matplotlib.pyplot as plt
import numpy as np
normaliser = np.count_nonzero(alignment_filtered, axis=0)
fig, ax = plt.subplots(figsize=(10,3))
ax.bar(range(phylop.shape[1]), np.nansum(phylop, axis=0)/phylop.shape[0], color = 'grey')
ax.margins(x=0, y=0)
#ax.set_ylim(3,6)
ax.set_xlim(450,550)
# %%
from ma_mapper import extract_bigwig
output=extract_bigwig.fetch_bigwig(metadata=coord_file, bigwig=bigwig_file, custom_id=True)
# %%
import numpy as np
import pyBigWig
def extract_bigwig(bigwig_file, chrom, start_list, end_list, strand):
    if isinstance(bigwig_file, str):
        bigwig = pyBigWig.open(bigwig_file)
    elif str(type(bigwig_file)) == "<class 'pyBigWig.bigWigFile'>":
        bigwig = bigwig_file

    bigwig_arrays = []
    for i in range(len(start_list)):
        # zero-based to one-based conversion [BED -> VCF]
        start = start_list[i]
        end = end_list[i]
        bigwig_array = bigwig.values(chrom, start, end, numpy = True)
        if strand == '-':
            bigwig_array = np.flip(bigwig_array)
        bigwig_arrays.append(bigwig_array)
    if strand == '-':
        bigwig_arrays.reverse()
    bigwig_out = np.concatenate(bigwig_arrays)
    return bigwig_out
# %%
metadata = pd.read_csv(coord_file, sep='\t')
metadata['meta_id'] = metadata.iloc[:,4]
meta_id = metadata.iloc[:,4].unique()

#%%
grouped = metadata.groupby('meta_id')

chrom_list = grouped.apply(lambda x: x.iloc[:,0].unique()[0], include_groups=False).tolist()
start_list = grouped.apply(lambda x: x.iloc[:,1].tolist(), include_groups=False).tolist()
end_list = grouped.apply(lambda x: x.iloc[:,2].tolist(), include_groups=False).tolist()
strand_list = grouped.apply(lambda x: x.iloc[:,3].unique()[0], include_groups=False).tolist()
# %%
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat 
with ProcessPoolExecutor(max_workers=40) as executor:
        results  = executor.map(extract_bigwig, repeat(bigwig_file), chrom_list, start_list,end_list,strand_list)
bigwig_out = []
for result in results:
    bigwig_out.append(result)
#%%
results = []
for i  in range(len(chrom_list)):
    print(i)
    results.append(extract_bigwig(bigwig_file, chrom_list[i], start_list[i], end_list[i], strand_list[i]))
# %%
import cProfile
#%%
def test():
    metadata = pd.read_csv(coord_file, sep='\t')
    metadata['meta_id'] = metadata.iloc[:,4]
    meta_id = metadata.iloc[:,4].unique()
    grouped = metadata.groupby('meta_id')

    chrom_list = grouped.apply(lambda x: x.iloc[:,0].unique()[0], include_groups=False).tolist()
    start_list = grouped.apply(lambda x: x.iloc[:,1].tolist(), include_groups=False).tolist()
    end_list = grouped.apply(lambda x: x.iloc[:,2].tolist(), include_groups=False).tolist()
    strand_list = grouped.apply(lambda x: x.iloc[:,3].unique()[0], include_groups=False).tolist()
    bigwig=pyBigWig.open(bigwig_file)
    results = []
    for i  in range(len(chrom_list)):
        results.append(extract_bigwig(bigwig, chrom_list[i], start_list[i], end_list[i], strand_list[i]))
#%%
cProfile.run('test()')
# %%
