#%%
from Bio.AlignIO import MafIO
import numpy as np
#%% setting
maf_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/241genomes/241-mammalian-2020v2b.maf.chr1'
target_species = 'Homo_sapiens'
chrom = 'chr1'
start = [4381114]
end = [4382190]
strand = '-'
coverage_count = False
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
                #test_list.append(np.array(list(seqrec.seq)))
                collector[seqrec.id.split('.')[0]]= np.array(list(seqrec.seq.lower()))
                print(seqrec.id,len(seqrec.seq.lower()),'\n',seqrec.seq.lower())
array_transposed=np.array(list(collector.values())).transpose()
count_arg = 'common'
output_array=[]
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
output_array
# %%
