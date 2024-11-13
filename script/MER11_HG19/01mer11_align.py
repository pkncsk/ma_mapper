#%% extract sequence
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import sequence_alignment
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
#%% load TE coords
hg19_coord_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11_HG19/t28_final_peaks_te_mer11_coords.bed'
hg19_coord_df=pd.read_csv(hg19_coord_file, sep='\t', header=None)
hg19_coord_df.columns=['chrom','start','end','subfamily','#peak']
hg19_coord_df['start'] = hg19_coord_df['start'] -1
#%% load repeatmasker table
repeatmasker_hg19 = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/repeatmasker/repeatmasker/hg19_repeatlib2014/hg19.fa.out.tsv'
repeatmasker_hg19_df=pd.read_csv(repeatmasker_hg19, sep='\t', index_col=0)
#%% make bed file for seqeuence extraction and alignment
prealign_metadata=hg19_coord_df.merge(repeatmasker_hg19_df, how='left',left_on=['chrom','start','end'], right_on=['genoName','genoStart','genoEnd'])[['chrom','start','end','subfamily','#peak','strand']]
prealign_metadata['name'] = 'MER11_'+prealign_metadata.index.astype(str)
prealign_metadata['score'] = 10
coord_file = prealign_metadata[['chrom','start','end','name','score','strand']]
#%% alignment
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from Bio.SeqRecord import SeqRecord
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import sequence_alignment

subfamily = 'MER11'
source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg19_fasta/hg19.fa'
metadata = coord_file
fasta_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11_HG19/mer11_w_peak.fasta'
#%% 
sequence_alignment.sequence_io(metadata,source_fasta,output_filepath =fasta_file, save_to_file= True,custom_id=True)
#%%
sequence_alignment.mafft_align(fasta_file, nthread = 6)
#%% parse alignment
import sys
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
subfamily = 'MER11'
from ma_mapper import mapper
alignment_file = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11_HG19/mer11_w_peak.fasta.aligned'
#%% parse and filter alignment into matrix
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file, col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
metadata_filtered
#%% merge prealignmetadata to add more information -> make plot annotation
metadata_w_info=metadata_filtered.merge(prealign_metadata[['subfamily','name','#peak']], how='left', on='name')

#%% make annotation
import numpy as np
subgroups = np.unique(metadata_w_info['subfamily'].astype(str))
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno=metadata_w_info['subfamily'].map(numerical_subgroup)
#%% plot
from ma_mapper import plots
import importlib
importlib.reload(plots)
plots.plot_experimental(
    show_alignment=True,
    alignment=alignment_filtered, 
    alignment_col='dna', 
    show_alignment_colbar=True,
    colorbar=True,  
    agg_major_tick=100, 
    annotation=True, 
    annotation_data=[subgroup_anno], 
    anno_col=[['red','yellow','blue']], 
    anno_title=['subfamily'],
    anno_cbar=True, 
    anno_cbar_label=[['MER11A','MER11B','MER11C']], )
# %% demonstration on how change to plot order need to be applied to every thing
metadata_sorted=metadata_w_info.sort_values('subfamily')
new_order=metadata_sorted.index
alignment_sorted=alignment_filtered[new_order]
subgroup_anno_sorted = subgroup_anno[new_order]

# %%
#%% plot
from ma_mapper import plots
import importlib
importlib.reload(plots)
plots.plot_experimental(
    show_alignment=True,
    alignment=alignment_sorted, 
    alignment_col='dna', 
    show_alignment_colbar=True,
    colorbar=True,  
    agg_major_tick=100, 
    annotation=True, 
    annotation_data=[subgroup_anno_sorted], 
    anno_col=[['red','yellow','blue']], 
    anno_title=['subfamily'],
    anno_cbar=True, 
    anno_cbar_label=[['MER11A','MER11B','MER11C']], )
# %% export tables
prealign_metadata.to_csv('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11_HG19/prealign_metadata.txt',index=False,sep='\t')
# just to prevent confusion, original_order refer to the order of alignment
metadata_w_info=metadata_w_info.drop(columns='original_order')
metadata_w_info.to_csv('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11_HG19/metadata.txt',index=False,sep='\t')

#%%
#export array
np.savetxt("/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/MER11_HG19/alignment_matrix_weakfilter.txt", alignment_filtered, delimiter="\t")
# %%
