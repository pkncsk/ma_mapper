#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import mapper
from ma_mapper import custom_cmap
#%%
subfamily = 'THE1C'
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
age_table = f'{config.te_age_folder}/{subfamily}.txt'
age_df = pd.read_csv(age_table, sep='\t')
internal_id_tbl = f'{config.internal_id_folder}/{subfamily}.txt'
internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
te_age_internal_id=internal_id_df.merge(age_df, on='internal_id',how='left')
#%%
age_default_id = pd.DataFrame()
age_default_id['internal_id'] = subfamily + '_' + te_age_internal_id.index.astype(str)
age_default_id['te_age'] = te_age_internal_id['te_age']
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered, age_table=age_default_id)
# %%
coord_file = f'{config.coord_internal_id_folder}/{subfamily}.txt'
maf_dir = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/241genomes/241-mammalian-2020v2b.maf'
e_value_table = f'{config.e_value_folder}/{subfamily}.txt'
#%%
internal_id_tbl = f'{config.internal_id_folder}/{subfamily}.txt'
internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
internal_id_df['meta_id'] = subfamily + '_' + internal_id_df.index.astype(str)
#%%
separated_maf = True
meta_id_list = metadata_age['name'].unique()
grouped = metadata_age.groupby('name', sort=False)
chrom_list = grouped.apply(lambda x: x.iloc[:,0].unique()[0], include_groups=False).tolist()
start_list = grouped.apply(lambda x: x.iloc[:,1].astype(int).tolist(), include_groups=False).tolist()
end_list = grouped.apply(lambda x: x.iloc[:,2].astype(int).tolist(), include_groups=False).tolist()
strand_list = grouped.apply(lambda x: x.iloc[:,5].unique()[0], include_groups=False).tolist()
maf_call_list = []
for chrom in chrom_list:
    if separated_maf == True:
        maf_file = f'{maf_dir}.{chrom}'
    else:
        maf_file = maf_dir
    maf_call_list.append(maf_file)
#%%
iupac_codes = {
        'A': {'A'},
        'C': {'C'},
        'G': {'G'},
        'T': {'T'},
        'R': {'A', 'G'},
        'Y': {'C', 'T'},
        'S': {'G', 'C'},
        'W': {'A', 'T'},
        'K': {'G', 'T'},
        'M': {'A', 'C'},
        'B': {'C', 'G', 'T'},
        'D': {'A', 'G', 'T'},
        'H': {'A', 'C', 'T'},
        'V': {'A', 'C', 'G'},
        'N': {'A', 'C', 'G', 'T'},
        '-': {'-'}
    }

def count_bases_with_ambiguity(sequences):
    base_counts = Counter()
    for sequence in sequences:
        for base in sequence:
            possible_bases = iupac_codes.get(base.upper(), {base.upper()})
            # Distribute counts among all possible bases
            for possible_base in possible_bases:
                base_counts[possible_base] += 1 / len(possible_bases)
    
    return base_counts
#%%
from Bio.AlignIO import MafIO
import numpy as np
target_species = 'Homo_sapiens'
import time
from ma_mapper import extract_maf
MafIO.MafIndex.get_spliced = extract_maf.get_spliced_mod
count_arg='base_count'
start_time = time.time()
test = []
species_list = ['Homo_sapiens', 'Gorilla_gorilla']
#for idx, meta_id in enumerate(meta_id_list):
meta_id = 'THE1C_100'
idx=100
print(idx, meta_id)
start=start_list[idx]
end=end_list[idx]
strand=strand_list[idx]
e_value_df = pd.read_csv(e_value_table, sep='\t')
maf_id = f'{target_species}.{chrom}'

index_maf = MafIO.MafIndex(f'{maf_file}.mafindex', maf_file, maf_id) 
n_strand = -1 if strand == '-' else 1
results =index_maf.get_spliced(start,end,n_strand)
    
#%%
if e_value_df is None and internal_id_df is None and species_list is None:
    collector = results
elif species_list is not None:
    pattern = '|'.join(species_list)  
    collector = results[results['seqid'].str.contains(pattern)]
else:
    if internal_id_df is not None:
        #subfamily = str.split(internal_id, sep ='_')[:-1]
        #internal_id_idx = str.split(internal_id, sep ='_')[-1]
        #_internal_id = internal_id_df[internal_id_df.index==int(internal_id_idx)]['internal_id'].values[0]
        _internal_id = internal_id_df[internal_id_df['meta_id'] == meta_id]['internal_id'].values[0]
    else:
        _internal_id = meta_id
    e_value_df['seqid'] = e_value_df.species + '.' + e_value_df.chr_code
    e_value_internal_id=e_value_df[e_value_df['internal_id'] == _internal_id]
    
    if e_value_internal_id.shape[0] <= 1:
        collector = results[results.seqid.str.contains('Homo_sapiens')]
    else:
        collector = results[results.seqid.isin(e_value_internal_id.seqid.values)]
#if count_arg == 'raw':
#    return collector
#return collector
#logger.info('pass1')
try:
    #logger.info(np.char.upper(results[results.seqid.str.contains('Homo_sapiens')]['seq'].to_list()))
    ref_alleles = np.char.upper(results[results.seqid.str.contains('Homo_sapiens')]['seq'].to_list())[0]
    #logger.info(ref_alleles)
except IndexError:
    print(f'IndexError:{meta_id}\t{maf_file}\t{maf_id}\t{chrom}\t{start}\t{end}\t{strand}')
    #print(collector)
    


array_transposed=np.array(collector['seq'].to_list()).transpose()
output_array=[]
for idx, pos_array in enumerate(array_transposed):
    frequencies =count_bases_with_ambiguity(np.char.upper(pos_array))
    ref_allele=ref_alleles[idx]
    #frequencies.pop('-', None) #->count deletion/insertion
    total = sum(frequencies.values())
    if count_arg == 'human_ref':
        alt_count = total - frequencies[ref_allele]
        alt_freq=alt_count/total
        if total == 1:
            alt_freq = np.nan
        output_array.append(alt_freq)
    elif count_arg in ['coverage', 'total_raw']:
        output_array.append(total)
    elif count_arg in ['common', 'common_raw']:
        common_allele = max(frequencies, key=frequencies.get)
        common_count = frequencies.get(common_allele)
        if count_arg == 'common':
            common_freq = common_count/total
            output_array.append(common_freq)
        else:
            output_array.append(common_count)
    elif count_arg == 'common_nogap':
        frequencies.pop('-', None)
        try:
            common_allele = max(frequencies, key=frequencies.get)
        except ValueError:
            print(f'ValueError:{meta_id}\t{maf_file}\t{maf_id}\t{chrom}\t{start}\t{end}\t{strand}')
        common_count = frequencies.get(common_allele)
        common_freq = common_count/total
        output_array.append(common_freq)
    elif count_arg == 'base_count':
        output_array.append(frequencies)
    elif count_arg in ['base_freq','base_freq_nogap']:
        if count_arg == 'base_freq_nogap':
            frequencies.pop('-', None)
        
        count_freq = {base: count / total for base, count in frequencies.items()}
        output_array.append(count_freq)
if isinstance(output_array, str):
    print('debug',meta_id,output_array)
time.time() - start_time
#%%
sequence_length = len(results.iloc[0]['seq'])
list_of_dfs = []
for i in range(sequence_length):
    # Extract the base at the i-th position for each row
    temp_df = results[['seqid']].copy()  # Start with the seqid column
    temp_df['seq'] = results['seq'].apply(lambda x: x[i])  # Add the base at the i-th position
    list_of_dfs.append(temp_df)
# %%
