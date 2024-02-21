#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/phd_project_development/dev/packaging_dir/ma_mapper/')
# %%
input_filepath = '/home/pc575/phd_project_development/data/repeatmasker/repeatmasker/hg38_repeatlib2014/hg38.fa.out.tsv'
rmskout_table=pd.read_csv(input_filepath, sep='\t', index_col = 0 )
# %%
main_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
subfamily = 'MER11A'
subfam_table = rmskout_table[rmskout_table.repName == subfamily]
subfam_table = subfam_table[subfam_table.genoName.isin(main_chr)]
subfam_coord = subfam_table[['genoName','genoStart','genoEnd','strand']]
subfam_coord['id'] = subfam_table.repName +'_'+ subfam_table.id.astype(str) + '_'+subfam_table.index.astype(str)
# %%
subfam_coord.to_csv('/home/pc575/phd_project_development/data/ma_mapper_output/mer11a_coord.txt', sep='\t', index= False)
#%% from age_div table
subfamily = 'MER11A'
input_filepath = '/home/pc575/phd_project_development/data/_ma_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/combined_age_div/combined_age_and_div.txt'
main_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
age_div_table = pd.read_csv(input_filepath, sep='\t')
subfam_table=age_div_table[age_div_table.repName == subfamily]
subfam_table = subfam_table[subfam_table.genoName.isin(main_chr)]
subfam_coord = subfam_table[['genoName','genoStart','genoEnd','strand']]
subfam_coord['id'] = 'n'+subfam_table.internal_id.astype(str)
#subfam_coord['id'] = 'n'+subfam_table.internal_id.astype(str) + '_' + subfam_table.index.astype(str)
# %%
subfam_coord.to_csv('/home/pc575/phd_project_development/data/_ma_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11a_coord_with_id.txt', sep='\t', index= False)
# %%
