#%%
from Bio.AlignIO import MafIO
import numpy as np
#%% setting
species_age_table_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/species241_info.tsv'
maf_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/241genomes/241-mammalian-2020v2b.maf.chr1'
target_species = 'Homo_sapiens'
chrom = 'chr1'
start = [4381114]
end = [4382190]
strand = '-'
coverage_count = False
import pandas as pd
species = pd.read_table(species_age_table_file, sep='\t')
age=50
age_arg = 'calibrate'
count_arg = 'common'
#%%

print(target_species,chrom, start, end, strand)

maf_id = f'{target_species}.{chrom}'

index_maf = MafIO.MafIndex(f'{maf_file}.mafindex', maf_file, maf_id) 
n_strand = -1 if strand == '-' else 1
results =index_maf.get_spliced(start,end,n_strand)
if age_arg is None:
    collector = {seqrec.id.split('.')[0]: np.array(list(seqrec.seq.lower())) for seqrec in results}
else:
    if age_arg == 'calibrate':
        collector = {seqrec.id.split('.')[0]: np.array(list(seqrec.seq.lower())) for seqrec in results if age >= species[species.meta_name == seqrec.id.split('.')[0]]['Estimated Time (MYA)'].to_list()[0]}
    elif age_arg == 'extract':
        for seqrec in results:
            species_age = species[species.meta_name == seqrec.id.split('.')[0]]['Estimated Time (MYA)'].to_list()[0]
            if species_age >= age:
                age = species_age



array_transposed=np.array(list(collector.values())).transpose()
output_array=[]
ref_alleles=collector['Homo_sapiens']
for idx,ref_pos in enumerate(array_transposed):
    (unique, counts) = np.unique(ref_pos, return_counts=True)
    frequencies = dict(zip(unique, counts))
    ref_allele=ref_alleles[idx]
    #frequencies.pop('-', None) #->count deletion/insertion
    total = sum(frequencies.values())
    if count_arg == 'human_ref':
        alt_count = total - frequencies[ref_allele]
        alt_freq=alt_count/total
        output_array.append(alt_freq)
    elif count_arg in ['coverage', 'total_raw']:
        output_array.append(total)
    elif count_arg in ['common', 'common_raw']:
        common_allele = max(frequencies, key=frequencies.get)
        common_count = frequencies[common_allele]
        if count_arg == 'common':
            common_freq = common_count/total
            output_array.append(common_freq)
        else:
            output_array.append(common_count)

# %%
output_array
# %%
collector = {seqrec.id.split('.')[0]: np.array(list(seqrec.seq.lower())) for seqrec in results}
#%%
subfamily = ['THE1C']
coord_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/coord_internal_id/'+subfamily[0]+'.txt'
metadata=pd.read_csv(coord_file, sep='\t')
maf_dir = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/241genomes/241-mammalian-2020v2b.maf'
# %%
results = []
for idx, row in metadata.iterrows():
    chrom = row.genoName
    start = [row.genoStart]
    end = [row.genoEnd]
    strand = row.strand
    print(target_species,chrom, start, end, strand)
    maf_id = f'{target_species}.{chrom}'
    maf_file = f'{maf_dir}.{chrom}'
    index_maf = MafIO.MafIndex(f'{maf_file}.mafindex', maf_file, maf_id) 
    n_strand = -1 if strand == '-' else 1
    results =index_maf.get_spliced(start,end,n_strand)
    if age_arg is None:
        collector = {seqrec.id.split('.')[0]: np.array(list(seqrec.seq.lower())) for seqrec in results}
    else:
        if age_arg == 'calibrate':
            collector = {seqrec.id.split('.')[0]: np.array(list(seqrec.seq.lower())) for seqrec in results if age >= species[species.meta_name == seqrec.id.split('.')[0]]['Estimated Time (MYA)'].to_list()[0]}
        elif age_arg == 'extract':
            for seqrec in results:
                species_age = species[species.meta_name == seqrec.id.split('.')[0]]['Estimated Time (MYA)'].to_list()[0]
                if species_age >= age:
                    age = species_age
    array_transposed=np.array(list(collector.values())).transpose()
    output_array=[]
    ref_alleles=collector['Homo_sapiens']
    for idx,ref_pos in enumerate(array_transposed):
        (unique, counts) = np.unique(ref_pos, return_counts=True)
        frequencies = dict(zip(unique, counts))
        ref_allele=ref_alleles[idx]
        #frequencies.pop('-', None) #->count deletion/insertion
        total = sum(frequencies.values())
        if count_arg == 'human_ref':
            alt_count = total - frequencies[ref_allele]
            alt_freq=alt_count/total
            output_array.append(alt_freq)
        elif count_arg in ['coverage', 'total_raw']:
            output_array.append(total)
        elif count_arg in ['common', 'common_raw']:
            common_allele = max(frequencies, key=frequencies.get)
            common_count = frequencies[common_allele]
            if count_arg == 'common':
                common_freq = common_count/total
                output_array.append(common_freq)
            else:
                output_array.append(common_count)
results.append(output_array)
# %%
