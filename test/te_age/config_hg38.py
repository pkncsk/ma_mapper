#%%
import pandas as pd
import numpy as np
# %% load rmskout 
input_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/repeatmasker/repeatmasker/hg38_repeatlib2014/hg38.fa.out.tsv'
#mouse
#rmskout_table=pd.read_csv(input_filepath, sep='\t', header=None,names = ['bin','swScore','milliDiv','milliDel','milliIns','genoName','genoStart','genoEnd','genoLeft','strand','repName','repClass','repFamily','repStart','repEnd','repLeft','id'])
#human
rmskout_table=pd.read_csv(input_filepath, sep='\t', index_col = 0, header = 0)
#filter table
main_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
discard_class = ['Simple_repeat', 'Satellite/telo','Low_complexity','snRNA','tRNA','Satellite','srpRNA', 'rRNA', 'scRNA', 'RNA', 'Satellite/centr', 'Satellite/acro',    'LTR/ERVL?','LTR/ERV1?','LTR?','LTR/Gypsy?','DNA/hAT?','RC?/Helitron?','DNA?/hAT-Tip100?','DNA?/PiggyBac?','SINE?/tRNA']
filtered_table=rmskout_table[(~rmskout_table['repClass'].isin(discard_class)) & (rmskout_table['genoName'].isin(main_chr))]

te_subfam_list=rmskout_table[(~rmskout_table['repClass'].isin(discard_class)) & (rmskout_table['genoName'].isin(main_chr))]['repName'].unique().tolist()
species = pd.read_table("/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/species241_info.tsv", sep='\t')
age_canon=np.sort(pd.unique(species['Estimated Time (MYA)']))
species['meta_name']=species['Species'].str.replace(' ','_')
species_table=species[species['Source']!='3. Zoonomia (not in alignment)']

age_ref_table_template = {'age': age_canon, 'representative': ['human', 'chimpanzee','gorilla','orangutan','gibbon','old-world monkey','new-world monkey','primates','flying lemur','treeshrew','rodents','modern eutheria','placental'],'representative_age': ['human (0)', 'chimpanzee (6.7)','gorilla (9.06)','orangutan (15.76)','gibbon (20.19)','old-world monkey (29.44)','new-world monkey (43.2)','primates (73.8)','flying lemur (76)','treeshrew (82)','rodents (90)','modern eutheria(96)','placental mammals (105+)'], 'note':['human', 'ape','ape','ape','gibbon','old-world monkey','new-world monkey','primates','flying lemur','primitive prosimian/primate','rodents', 'Boreoeutheria', 'Atlantogenata']}
#%%
subfam_tally=filtered_table.groupby('repName').count()['swScore']
subfam1k=subfam_tally[subfam_tally<1000].index
subfam5k=subfam_tally[(subfam_tally<=5000) & (subfam_tally>1000)].index
subfam10k=subfam_tally[(subfam_tally<=10000) & (subfam_tally>5000)].index
subfam50k=subfam_tally[(subfam_tally<=50000) & (subfam_tally>10000)].index
subfam100k=subfam_tally[(subfam_tally<=100000) & (subfam_tally>50000)].index
subfam_large =subfam_tally[subfam_tally>100000].index
#['AluJb', 'AluJr', 'AluSx', 'AluSx1', 'AluSz', 'L2a', 'L2c', 'MIR','MIRb', 'MIRc'],
#%% app setting CURRENT
target_species = 'Homo_sapiens'
#subfamily_list = te_subfam_list
#subfamily_list = te_subfam_list[0:50]
#subfamily_list = te_subfam_list[50:100]
#subfamily_list = te_subfam_list[200:400]
#subfamily_list = te_subfam_list[400:600]
#subfamily_list = te_subfam_list[600:800]
#subfamily_list = te_subfam_list[800:1000]
#subfamily_list = te_subfam_list[100:]
#subfamily_list = te_subfam_list
#subfamily_list =subfam_tally.index
#subfamily_list = subfam100k
subfamily_list =['MIRc']
main_folder = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker/' 
internal_id_folder = main_folder + '/internal_id'
coord_internal_id_folder = f'{main_folder}/old_result_redo/coord_internal_id'
e_value_folder = f'{main_folder}/e_value'
te_age_folder = main_folder + '/age'
te_age_bed_folder = f'{main_folder}/age_evi_filtered_bed'
te_age_human_insertion_folder = main_folder + '/age/human_insertion'
te_age_segmental_folder = main_folder + '/age/segmental'
te_age_unclassified = main_folder + '/age/unclassified'
te_age_dist_fig_folder = main_folder + '/age_dist_fig'
te_length_folder = main_folder + '/length'
te_tag_folder = main_folder + '/tags'
outlier_surrounding_folder = main_folder + '/outlier_surrounding_age'
te_seq_folder = main_folder + '/seq'
te_seq_rmsk_folder = main_folder + '/seq_rmsk'
te_alignment_folder = main_folder + '/alignment'
te_alignment_rmsk_folder = main_folder + '/alignment_rmsk'
parsed_alignment_folder = main_folder + '/parsed_alignment'
parsed_alignment_filter_folder = main_folder + '/parsed_alignment_filter'
filtered_alignment_folder = main_folder + '/filtered_alignment'
te_alignment_fig_folder = main_folder + '/alignment_fig'
kimura_distance_folder = main_folder + '/div_redo'
combined_age_div_folder = main_folder + '/combined_age_div'
#subfamily_list = [['MER11A','MER11B','MER11C']]
#subfamily_list = ['MER11A','MER11B','MER11C']
#%%
subfam_tally=filtered_table.groupby('repName').count()['swScore']
# %%
