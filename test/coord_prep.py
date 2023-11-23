#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/phd_project_development/dev/packaging_dir/ma_mapper/')
# %%
input_filepath = '/home/pc575/phd_project_development/data/repeatmasker/repeatmasker/hg38_repeatlib2014/hg38.fa.out.tsv'
rmskout_table=pd.read_csv(input_filepath, sep='\t', index_col = 0 )
# %%
subfamily = 'MER11A'
subfam_table = rmskout_table[rmskout_table.repName == subfamily]
subfam_coord = subfam_table[['genoName','genoStart','genoEnd','strand']]
subfam_coord['id'] = subfam_table.repName +'_'+ subfam_table.id.astype(str) + '_'+subfam_table.index.astype(str)
# %%
subfam_coord.to_csv('/home/pc575/phd_project_development/data/ma_mapper_output/mer11a_coord.txt', sep='\t', index= False)
#%%