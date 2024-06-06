#%% extract sequence
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
subfamily = ['pcg_sliced']
from ma_mapper import sequence_alignment
coord_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/coord_internal_id/'+subfamily[0]+'.txt'
#%% alignment
import sys
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import sequence_alignment
subfamily = ['pcg_sliced']
source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_fasta/hg38.fa'
metadata = coord_file
fasta_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/alignment/'+subfamily[0]+'.fasta'
#%%
sequence_alignment.sequence_io(metadata,source_fasta,output_filepath =fasta_file, save_to_file= True,custom_id= True)
#%%
#sequence_alignment.mafft_align(fasta_file, nthread = 40)
#%% parse alignment
import sys
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
subfamily = ['pcg_sliced']
from ma_mapper import mapper
alignment_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/alignment/'+subfamily[0]+'.fasta'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered)
#%%
from ma_mapper import plot
plot.overlay_plot([alignment_filtered], metadata= metadata_age,age_annotation=False)
#%%
alignment_extend, metadata_filtered= mapper.parse_and_filter(alignment_file, extension_length=500)
#%%
plot.overlay_plot([alignment_extend], metadata= metadata_age,age_annotation=False)
#%%
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
records = (r for r in SeqIO.parse(alignment_file, "fasta"))
# %%
for i, value in enumerate(records):
    if len(value)!=1000:
        print(i, len(value),value)
# %%
