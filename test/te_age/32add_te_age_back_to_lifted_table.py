#%%
import pandas as pd
#%%
combined_te_age_lifted=pd.read_csv('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/combined_te_age_hg19_evi.bed',header=None, sep='\t')
combined_te_age_lifted.columns = ['chrom','start','end','name','del','strand']
# %%
combined_te_age = pd.read_csv('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/age_evi_filtered_bed/combined_te_age.txt',header=None, sep='\t')
combined_te_age.columns = ['chrom','start','end','name','te_age','strand']
# %%
liftover_corrected=combined_te_age_lifted.merge(combined_te_age, how='left', on=['name']).drop('del', axis=1)
liftover_corrected=liftover_corrected[['chrom_x','start_x','end_x','name','te_age','strand_x']]
# %%
liftover_corrected.columns = ['chrom','start','end','name','te_age','strand']
#liftover_corrected_bed =liftover_corrected[['chrom','start','end','name','te_age','strand']]
# %%
unique_te_age=combined_te_age[['name', 'te_age']].drop_duplicates()
# %%
combined_hg19=combined_te_age_lifted.merge(unique_te_age, how='left',on='name')
# %%
combined_hg19_bed = combined_hg19[['chrom','start','end','name','te_age','strand']]
# %%
combined_hg19_bed.to_csv('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/age_evi_filtered_bed/combined_te_age_hg19_evi.txt', sep='\t', index=False, header=False)
# %%
