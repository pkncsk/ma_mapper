#%%
import pandas as pd
import numpy as np
#%%
def repeatmasker_prep(repeatmasker_filepath):
    rmskout_table=pd.read_csv(repeatmasker_filepath, sep='\t', index_col = 0, header = 0)
    #filter table
    main_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
    discard_class = ['Simple_repeat', 'Satellite/telo','Low_complexity','snRNA','tRNA','Satellite','srpRNA', 'rRNA', 'scRNA', 'RNA', 'Satellite/centr', 'Satellite/acro',    'LTR/ERVL?','LTR/ERV1?','LTR?','LTR/Gypsy?','DNA/hAT?','RC?/Helitron?','DNA?/hAT-Tip100?','DNA?/PiggyBac?','SINE?/tRNA']
    filtered_table=rmskout_table[(~rmskout_table['repClass'].isin(discard_class)) & (rmskout_table['genoName'].isin(main_chr))]
    return filtered_table

def divergence_table_prep(divergence_table_filepath):
    species = pd.read_table(divergence_table_filepath, sep='\t')
    species['meta_name']=species['Species'].str.replace(' ','_')
    species_table=species[species['Source']!='3. Zoonomia (not in alignment)']
    return species

def age_canon():
    age_canon = [0,6.7,9.06,15.76,20.19,29.44,43.2,73.8,76,82,90,96,105,np.nan]
    return age_canon

def age_reference_table():
    age_canon = [0,6.7,9.06,15.76,20.19,29.44,43.2,73.8,76,82,90,96,105,np.nan]
    age_ref_table_template = {'age': age_canon, 'representative':
                          ['human', 
                           'chimpanzee',
                           'gorilla',
                           'orangutan',
                           'gibbon',
                           'old-world monkey',
                           'new-world monkey',
                           'primates',
                           'close primate relatives',
                           'early primate ancestor',
                           'rodents',
                           'modern eutheria',
                           'placental',
                           'undefined'],
                           'representative_age': 
                           ['human (0)', 
                            'chimpanzee (6.7)',
                            'gorilla (9.06)',
                            'orangutan (15.76)',
                            'gibbon (20.19)',
                            'old-world monkey (29.44)',
                            'new-world monkey (43.2)',
                            'primates (73.8)',
                            'close primate relatives (76)',
                            'early primate ancestor (82)',
                            'rodents (90)',
                            'modern eutheria (96)',
                            'placental mammals (105+)',
                            'Undefined (NA)'], 
                            'note':
                            ['human', 
                             'ape',
                             'ape',
                             'ape',
                             'gibbon',
                             'old-world monkey',
                             'new-world monkey',
                             'primates',
                             'colugos (flying Lemurs)',
                             'primitive prosimian/primate',
                             'rodents', 
                             'Boreoeutheria', 
                             'Atlantogenata',
                             'Undefined (NA)']}
    return age_ref_table_template
#%%
def age_reference_table_reorg():
    age_canon_reorg = [0,6.7,9.06,15.76,20.19,29.44,43.2,73.8,76,96,105,np.nan]
    age_ref_table_reorg = {'age': age_canon_reorg, 'representative': ['human', 'chimpanzee','gorilla','orangutan','gibbon','old-world monkey','new-world monkey','primates','close primate relatives','modern eutheria','placental','undefined'],'representative_age': ['human (0)', 'chimpanzee (6.7)','gorilla (9.06)','orangutan (15.76)','gibbon (20.19)','old-world monkey (29.44)','new-world monkey (43.2)','primates (73.8)','close primate relatives (76)','modern eutheria (96)','placental mammals (105+)','undefined (NA)'], 'note':['human', 'ape','ape','ape','gibbon','old-world monkey','new-world monkey','primates','colugos (flying Lemurs)/primitive prosimian/primate','rodents/Boreoeutheria', 'Atlantogenata','Undefinied']}
    return age_ref_table_reorg