#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import extract_bigwig
from ma_mapper import extract_bam
from ma_mapper import mapper
#%%
subfamily = ['THE1C']
alignment_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/alignment/'+subfamily[0]+'.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered)
# %%
coord_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/coord_internal_id/'+subfamily[0]+'.txt'
bigwig_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_bigwig/241-mammalian-2020v2.bigWig'
#%%
phylop=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig')

# %%
from ma_mapper import plot
plot.overlay_plot([alignment_filtered,phylop], metadata= metadata_age,data_vmin=[-1.0],data_vmax=[1.0], nucleotide_color='white',data_cmap=['RdBu'], plot_title='THE1C:Zoonomia_phyloP', show_data_legend=True,data_label=['phyloP'])
#%%
from ma_mapper import extract_bigwig

# %%
import numpy as np
def extract_bigwig(bigwig_file, chrom, start_list, end_list, strand):
    bigwig_arrays = []
    for i in range(len(start_list)):
        # zero-based to one-based conversion [BED -> VCF]
        start = start_list[i]
        end = end_list[i]

        bigwig = pyBigWig.open(bigwig_file)
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
import pyBigWig   
chrom_list = []
start_list = []
end_list = []
strand_list = []    
for uniq_meta_id in meta_id:
        metadata_by_id = metadata[metadata.meta_id == uniq_meta_id]
        chrom = metadata_by_id.iloc[:,0].unique()[0]
        chrom_list.append(chrom)
        start_list.append(metadata_by_id.iloc[:,1].to_list())
        end_list.append(metadata_by_id.iloc[:,2].to_list())
        strand_list.append(metadata_by_id.iloc[:,3].unique()[0])
#%%
extract_bigwig(bigwig_file, chrom_list[0], start_list[0], end_list[0], strand_list[0])
# %%
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat 
with ProcessPoolExecutor(max_workers=40) as executor:
        results  = executor.map(extract_bigwig, repeat(bigwig_file), chrom_list, start_list,end_list,strand_list)
bigwig_out = []
for result in results:
    bigwig_out.append(result)
# %%
