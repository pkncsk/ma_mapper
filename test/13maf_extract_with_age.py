#%%
from Bio.AlignIO import MafIO
import numpy as np
import pandas as pd
#%%
#%% from age_div table
subfamily = ['MER11A','MER11B','MER11C']
input_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/combined_age_div/combined_age_and_div.txt'
main_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
age_div_table = pd.read_csv(input_filepath, sep='\t')
subfam_table=age_div_table[age_div_table.repName.isin(subfamily)]
subfam_table = subfam_table[subfam_table.genoName.isin(main_chr)]
subfam_table['id'] = subfam_table.repName+'_'+subfam_table.internal_id.astype(str)
#%%
species = pd.read_table("/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/species241_info.tsv", sep='\t')
#%% setting
maf_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/241genomes/241-mammalian-2020v2b.maf.chr1'
target_species = 'Homo_sapiens'
chrom = 'chr1'
start = [4380614]
end = [4381114]
strand = '-'
coverage_count = False
age = 96.00
#%%
print(target_species,chrom, start, end, strand)
maf_id = target_species+'.'+chrom
from Bio.AlignIO import MafIO
index_maf = MafIO.MafIndex(maf_file+".mafindex", maf_file, maf_id) 
if strand =='-':
    n_strand = -1
    results =index_maf.get_spliced(start,end,n_strand)
else:
    n_strand = 1
    results =index_maf.get_spliced(start,end,n_strand)
collector = {}
for seqrec in results:
    #print(seqrec.id.split('.')[0])
    if seqrec.id.split('.')[0] == target_species:  
        collector[seqrec.id.split('.')[0]]= np.array(list(seqrec.seq.lower()))
        print(seqrec.id.split('.')[0],len(seqrec.seq.lower()),'\n',seqrec.seq.lower())
        for seqrec in results:
            if seqrec.id.split('.')[0] != target_species:
                if age is not None:
                    species_age = species[species.meta_name == seqrec.id.split('.')[0]]['Estimated Time (MYA)'].to_list()[0]
                    print(seqrec.id.split('.')[0])
                    if age >= species_age:
                        #test_list.append(np.array(list(seqrec.seq)))
                        collector[seqrec.id.split('.')[0]]= np.array(list(seqrec.seq.lower()))
                        print(seqrec.id.split('.')[0],species_age,len(seqrec.seq.lower()),'\n',seqrec.seq.lower())
                    else:
                        continue
                else:
                    collector[seqrec.id.split('.')[0]]= np.array(list(seqrec.seq.lower()))
                    print(seqrec.id.split('.')[0],len(seqrec.seq.lower()),'\n',seqrec.seq.lower())
array_transposed=np.array(list(collector.values())).transpose()
alt_freq_array=[]
for ref_pos in array_transposed:
    (unique, counts) = np.unique(ref_pos, return_counts=True)
    frequencies = dict(zip(unique, counts))
    ref_allele=ref_pos[0]
    total = 0
    alt_count  = 0
    for key in frequencies:
        if key != '-':
            total = total+frequencies[key]
            if key != ref_allele:
                alt_count =alt_count+frequencies[key]
    if coverage_count is True:
        alt_freq = total
    else:
        alt_freq=alt_count/total
    alt_freq_array.append(alt_freq)
# %%
count_arg = 'human_ref'
output_array  = []
for ref_pos in array_transposed:
    (unique, counts) = np.unique(ref_pos, return_counts=True)
    frequencies = dict(zip(unique, counts))
    ref_allele=ref_pos[0]
    frequencies.pop('-', None)
    total = sum(frequencies.values())
    if count_arg == 'human_ref':
        alt_count = total - frequencies[ref_allele]
        alt_freq=alt_count/total
        output_array.append(alt_freq)
    elif count_arg == 'coverage':
        output_array.append(total)
    elif count_arg == 'common':
        common_allele = max(frequencies, key=frequencies.get)
        common_count = frequencies[common_allele]
        common_freq = common_count/total
        output_array.append(common_freq)
# %%
for i in output_array:
    print(i)
# %%
