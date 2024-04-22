#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import fetch_sequence
from ma_mapper import mafft_align
import argparse
#%%
#parser = argparse.ArgumentParser(prog='alignment pipeline',description='take subfamily input and run the pipe')
#parser.add_argument('subfamily', type = str)
#args = parser.parse_args()

#%%
rmsk_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/repeatmasker/repeatmasker/hg38_repeatlib2014/hg38.fa.out.tsv'
rmskout_table=pd.read_csv(rmsk_filepath, sep='\t', index_col = 0 )
discard_class = ['Simple_repeat','Satellite', 'rRNA', 'scRNA','srpRNA','tRNA','snRNA','Low_complexity',] 
main_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
#%% 
te_subfam_list=rmskout_table[(~rmskout_table['repClass'].isin(discard_class)) & (rmskout_table['genoName'].isin(main_chr))]['repName'].unique().tolist()
#NOTE: subfamily list just in case
#%% from age_div table
#subfamily = args.subfamily #NOTE: subfamily name here
subfamily = 'THE1C'
input_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/combined_age_div/combined_age_and_div.txt'
age_div_table = pd.read_csv(input_filepath, sep='\t')
subfam_table=age_div_table[age_div_table.repName == subfamily]
subfam_table = subfam_table[subfam_table.genoName.isin(main_chr)]
subfam_coord = subfam_table[['genoName','genoStart','genoEnd','strand']]
subfam_coord['id'] = 'n'+subfam_table.internal_id.astype(str)
#%%
coord_folder = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo' #NOTE: a folder to store coordinate files
coord_filepath = coord_folder+'/'+subfamily+'.txt'
subfam_coord.to_csv(coord_filepath, sep='\t', index= False)
#%% MULTIPLE ALIGNMENT
source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_fasta/hg38.fa'
metadata = coord_filepath
fasta_folder = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo' #NOTE: a folder to store fasta files and alignments
fasta_file = fasta_folder +'/'+subfamily+'.fasta' 
fetch_sequence.fetch_sequence(metadata,source_fasta,output_filepath =fasta_file, custom_id= True, save_to_file=True)

#some parameters can be changed from here
# nthread, nthreadtb, nthreadit
# additional commands can be inserted to mafft_arg = '(insert parameter here)'
# output_filepath can be set to store alignment files separately from fasta
mafft_align.mafft_wrapper(fasta_file, nthread = 40)
#%%