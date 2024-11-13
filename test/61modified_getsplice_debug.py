#%%
from Bio.AlignIO import MafIO
from itertools import repeat 
import os
import sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd

#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import mapper
from ma_mapper import custom_cmap
from ma_mapper import sequence_alignment
#%%
subfamily = 'MER11A'
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
age_table = f'{config.te_age_folder}/{subfamily}.txt'
age_df = pd.read_csv(age_table, sep='\t')
#%%
internal_id_tbl = f'{config.internal_id_folder}/{subfamily}.txt'
internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
internal_id_sort = internal_id_df.sort_values('rmsk_index')

te_age_internal_id=internal_id_sort.merge(age_df, on='internal_id', how='left')
internal_id_sort['meta_id'] =subfamily + '_' + internal_id_sort.index.astype(str)
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
#%%
separated_maf = True
meta_id_list = metadata_age['name'].unique()
grouped = metadata_age.groupby('name', sort=False)
chrom_list = grouped.apply(lambda x: x.iloc[:,0].unique()[0], include_groups=False).tolist()
start_list = grouped.apply(lambda x: x.iloc[:,1].astype(int).tolist(), include_groups=False).tolist()
end_list = grouped.apply(lambda x: x.iloc[:,2].astype(int).tolist(), include_groups=False).tolist()
strand_list = grouped.apply(lambda x: x.iloc[:,4].unique()[0], include_groups=False).tolist()
maf_call_list = []
for chrom in chrom_list:
    if separated_maf == True:
        maf_file = f'{maf_dir}.{chrom}'
    else:
        maf_file = maf_dir
    maf_call_list.append(maf_file)
#%%
from Bio.AlignIO import MafIO
import numpy as np
target_species = 'Homo_sapiens'
import time
count_arg='base_count'
start_time = time.time()
test = []
#%%
def get_spliced_mod(self, starts, ends, strand=1):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from collections import Counter, defaultdict

    # Dictionary for IUPAC ambiguity codes for 2-base combinations
    iupac_code = {
        frozenset(['A', 'G']): 'R', frozenset(['C', 'T']): 'Y',
        frozenset(['G', 'C']): 'S', frozenset(['A', 'T']): 'W',
        frozenset(['G', 'T']): 'K', frozenset(['A', 'C']): 'M',
        frozenset(['C', 'G', 'T']): 'B', frozenset(['A', 'G', 'T']): 'D',
        frozenset(['A', 'C', 'T']): 'H', frozenset(['A', 'C', 'G']): 'V',
        frozenset(['A', 'C', 'G', 'T']): 'N'
    }

    def convert_to_iupac(sequence):
        unique_bases = frozenset(sequence)
        if len(unique_bases) == 1:
            return sequence[0].upper()  
        return iupac_code.get(unique_bases, 'N')  # Default to 'N' for any unhandled cases

    def process_sequences_work(data):
        for species, positions in data.items():
            for pos, sequence in positions.items():
                sequence = sequence.upper()  # Convert the sequence to uppercase
                if len(sequence) > 1:
                    # Count the occurrences of each base
                    base_counts = Counter(sequence)
                    most_common_bases = base_counts.most_common()
                    max_count = most_common_bases[0][1]
                    consensus_bases = [base for base, count in most_common_bases if count == max_count]
                    # Convert the most common bases to IUPAC code
                    new_base = convert_to_iupac(consensus_bases)
                    
                    # Update the sequence with the new base
                    data[species][pos] = new_base
                else:
                    # Keep the single base as it is
                    data[species][pos] = sequence.upper()
        return data
    def process_sequence_localized(sequence):
        sequence = sequence.upper()
        base_counts = Counter(sequence)
        most_common_bases = base_counts.most_common()
        max_count = most_common_bases[0][1]
        consensus_bases = [base for base, count in most_common_bases if count == max_count]
        new_base = convert_to_iupac(consensus_bases)
        return new_base

    if strand not in (1, -1):
        raise ValueError("Strand must be 1 or -1, got %s" % strand)
    fetched = list(self.search(starts, ends))
    #return fetched
    expected_letters = sum(end - start for start, end in zip(starts, ends))
    if len(fetched) == 0:
        return [SeqRecord(Seq("N" * expected_letters), id=self._target_seqname)]
    all_seqnames = {sequence.id for multiseq in fetched for sequence in multiseq}
    split_by_position = {seq_name: {} for seq_name in all_seqnames}
    split_by_position
    total_rec_length = 0
    ref_first_strand = None
    #return all_seqnames
    first_step = []
    second_step = []
    third_step = []
    for multiseq in fetched:
        start_time = time.time()
        for seqrec in multiseq:
            #if seqrec.id == self._target_seqname:
            if seqrec.id == self._target_seqname:
                try:
                    if ref_first_strand is None:
                        ref_first_strand = seqrec.annotations["strand"]

                        if ref_first_strand not in (1, -1):
                            raise ValueError("Strand must be 1 or -1")
                    elif ref_first_strand != seqrec.annotations["strand"]:
                        raise ValueError(
                            "Encountered strand='%s' on target seqname, "
                            "expected '%s'"
                            % (seqrec.annotations["strand"], ref_first_strand)
                        )
                except KeyError:
                    raise ValueError(
                        "No strand information for target seqname (%s)"
                        % self._target_seqname
                    ) from None

                rec_length = len(seqrec)
                rec_start = seqrec.annotations["start"]
                ungapped_length = seqrec.annotations["size"]
                rec_end = rec_start + ungapped_length - 1
                total_rec_length += ungapped_length
                
                for seqrec in multiseq:
                    for pos in range(rec_start, rec_end + 1):
                        split_by_position[seqrec.id][pos] = ""

                break 
            else:
                raise ValueError(
                    "Did not find %s in alignment bundle" % (self._target_seqname,)
                )
        first_step.append(time.time()-start_time)
        start_time = time.time()
        real_pos = rec_start
        
        edit_id = []
        edit_pos = []
        for gapped_pos in range(rec_length):
            previous_id = ''
            for seqrec in multiseq:
                
                if seqrec.id == self._target_seqname:
                    track_val = seqrec.seq[gapped_pos]
                
                
                split_by_position[seqrec.id][real_pos] += seqrec.seq[gapped_pos]
                if previous_id == seqrec.id:
                        edit_id.append(seqrec.id)
                        edit_pos.append(real_pos)
                previous_id = seqrec.id
            if track_val != "-" and real_pos < rec_end:
                real_pos += 1
        # Debugging: Print lengths of sequences in split_by_position
        second_step.append(time.time()-start_time)
        start_time = time.time()
        for i in range(len(edit_id)):
            #print(f"{edit_id[i]} position {edit_pos[i]}")
            _sequence=split_by_position[edit_id[i]][edit_pos[i]]
            new_sequence=process_sequence_localized(_sequence)
            split_by_position[edit_id[i]][edit_pos[i]] = new_sequence


        
        
        #split_by_position = process_sequences_work(split_by_position) 
        third_step.append(time.time()-start_time)
        #return(split_by_position)
        if len(split_by_position[self._target_seqname]) != total_rec_length:
            raise ValueError(
                "Target seqname (%s) has %s records, expected %s"
                % (
                    self._target_seqname,
                    len(split_by_position[self._target_seqname]),
                    total_rec_length,
                )
            )
    #return split_by_position
    start_time = time.time()
    realpos_to_len = {
        pos: len(gapped_fragment)
        for pos, gapped_fragment in split_by_position[self._target_seqname].items()
        if len(gapped_fragment) > 1
    }
    print('extract data:', sum(first_step))
    print('position filling:', sum(second_step))
    print('remove dupes and calculate consensus:', sum(third_step))
    #
    # splice together the exons
    #start_time = time.time()
    seqid_list = []
    seq_list = []
    
    for seqid in all_seqnames:
        seq_split = split_by_position[seqid]
        seq_splice = []
        filler_char = "N" if seqid == self._target_seqname else "-"
        append = seq_splice.append
        for exonstart, exonend in zip(starts, ends):
            for real_pos in range(exonstart, exonend):
                if real_pos in seq_split:
                    append(seq_split[real_pos])
                elif real_pos in realpos_to_len:
                    append(filler_char * realpos_to_len[real_pos])
                else:
                    append(filler_char)

        seqid_list.append(seqid)
        seq_list.append(Seq("".join(seq_splice))) 
    print('splicing:', time.time()-start_time)
    start_time = time.time()
    if len(seq_list[seqid_list.index(self._target_seqname)].replace("-", "")) != expected_letters:
        raise ValueError(
            "Returning %s letters for target seqname (%s), expected %s"
            % (
                len(seq_list[seqid_list.index(self._target_seqname)].replace("-", "")),
                self._target_seqname,
                expected_letters,
            )
        )
    
    ref_subseq_len = len(seq_list[seqid_list.index(self._target_seqname)])
    for seqid, seq in zip(seqid_list, seq_list):
        if len(seq) != ref_subseq_len:
            raise ValueError(
                "Returning length %s for %s, expected %s"
                % (len(seq), seqid, ref_subseq_len)
            )
    
    # Create a DataFrame
    df = pd.DataFrame({
        'seqid': seqid_list,
        'seq': [seq.reverse_complement() if strand != ref_first_strand else seq for seq in seq_list]
    })
    print('dataframe creation:', time.time()-start_time)
    return df
#%%
MafIO.MafIndex.get_spliced = get_spliced_mod
#%%
for idx, meta_id in enumerate(['MER11A_1']):
        start_time = time.time()
        start=start_list[idx]
        end=end_list[idx]
        strand=strand_list[idx]
        chrom = chrom_list[idx]
        #e_value_df = pd.read_csv(e_value_table, sep='\t')
        maf_id = f'{target_species}.{chrom}'
        maf_file = maf_call_list[idx]
        start_flanked=[min(start)-5000] + start + [max(end)]
        end_flanked = [min(start)] + end + [max(end)+5000]
        index_maf = MafIO.MafIndex(f'{maf_file}.mafindex', maf_file, maf_id) 
        n_strand = -1 if strand == '-' else 1
        results =index_maf.get_spliced(start,end,n_strand)
        print(idx, meta_id, start_flanked, end_flanked, strand, time.time()-start_time)
        break
# %%
e_value_df = pd.read_csv(e_value_table, sep= '\t')
e_value_df['seqid'] = e_value_df.species + '.' + e_value_df.chr_code
# %%
# %%
internal_id = 'MER11A_1'
subfamily = str.split(internal_id, sep ='_')[:-1]
internal_id_idx = str.split(internal_id, sep ='_')[-1]
_internal_id = internal_id_sort[internal_id_sort.meta_id == internal_id]['internal_id'].values[0]
e_value_internal_id=e_value_df[e_value_df.internal_id == _internal_id]
# %%
results[results.seqid.isin(e_value_internal_id.seqid.values)]
# %%
results[results.seqid.str.contains('Homo_sapiens')]
# %%
from collections import Counter

# Define IUPAC ambiguity codes and their corresponding sets
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
    'N': {'A', 'C', 'G', 'T'}
}

def count_bases_with_ambiguity(sequences):
    # Initialize a Counter for base frequencies
    base_counts = Counter()

    for sequence in sequences:
        for base in sequence:
            possible_bases = iupac_codes.get(base.upper(), {base.upper()})
            # Distribute counts among all possible bases
            for possible_base in possible_bases:
                base_counts[possible_base] += 1 / len(possible_bases)
    
    return base_counts

# Example sequences
sequences = [
    'AGCT',
    'AGGT',
    'AGTT',
    'AGAT',
    'R',     # Represents A or G
    'Y',     # Represents C or T
    'W'      # Represents A or T
]

# Count base frequencies
base_frequencies = count_bases_with_ambiguity(sequences)

# %%
results[results.seqid.str.contains('Homo_sapiens')]
# %%
test=np.char.upper(results[results.seqid.str.contains('Homo_sapiens')]['seq'].to_list())[0]
# %%
count_bases_with_ambiguity(test)
# %%
iupac_code = {
    frozenset(['A', 'G']): 'R', frozenset(['C', 'T']): 'Y',
    frozenset(['G', 'C']): 'S', frozenset(['A', 'T']): 'W',
    frozenset(['G', 'T']): 'K', frozenset(['A', 'C']): 'M',
    frozenset(['C', 'G', 'T']): 'B', frozenset(['A', 'G', 'T']): 'D',
    frozenset(['A', 'C', 'T']): 'H', frozenset(['A', 'C', 'G']): 'V',
    frozenset(['A', 'C', 'G', 'T']): 'N'
}

def convert_to_iupac(sequence):
    unique_bases = frozenset(sequence)
    if len(unique_bases) == 1:
        return sequence[0].upper()  
    return iupac_code.get(unique_bases, 'N')  # Default to 'N' for any unhandled cases

def process_sequence_localized(sequence):
    sequence = sequence.upper()
    filtered_sequence = [base for base in sequence if base != '-']

    #base_counts = Counter(sequence)
    #most_common_bases = base_counts.most_common()
    #max_count = most_common_bases[0][1]
    #consensus_bases = [base for base, count in most_common_bases if count == max_count]
    new_base = convert_to_iupac(filtered_sequence)
    return new_base
# %%
