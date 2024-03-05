#%%
from Bio.AlignIO import MafIO
import numpy as np
#%% setting
maf_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/241genomes/241-mammalian-2020v2b.maf.chrY'
target_species = 'Homo_sapiens'
chrom = 'chrY'
start = [57150414]
end = [57151452]
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
