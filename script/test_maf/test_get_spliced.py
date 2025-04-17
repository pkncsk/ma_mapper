#%%
from math import exp
import sys
sys.path.append('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/ma_mapper/')
from ma_mapper import extract_maf
from ma_mapper import mapper
import pandas as pd
import os
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat 
#%%
name = 'test'
chrom = 'chr1'#'chr2'
maf_file = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447//{chrom}.maf'

start = [119563] #[156425872]
end = [119944] #[156425923]
strand = '-'
target_species = 'hg38'

#%%
from Bio.AlignIO import MafIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
def get_spliced_mod(self, starts, ends, strand=1):
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
    
    def process_sequence_localized(sequence):
        sequence = sequence.upper()
        filtered_sequence = [base for base in sequence if base != '-']

        #base_counts = Counter(sequence)
        #most_common_bases = base_counts.most_common()
        #max_count = most_common_bases[0][1]
        #consensus_bases = [base for base, count in most_common_bases if count == max_count]
        new_base = convert_to_iupac(filtered_sequence)
        return new_base

    if strand not in (1, -1): 
        raise ValueError("Strand must be 1 or -1, got %s" % strand)
    fetched = list(self.search(starts, ends))
    return fetched
    
MafIO.MafIndex.get_spliced = get_spliced_mod
#%%
from collections import Counter, defaultdict
maf_id = f'{target_species}.{chrom}'
mafindex_filedir = '.'.join(str.split(maf_file, sep='.')[:-1])
mafindex_filepath = f'{mafindex_filedir}.mafindex'
index_maf = MafIO.MafIndex(mafindex_filepath, maf_file, maf_id) 
n_strand = -1 if strand == '-' else 1
fetched =index_maf.get_spliced(start,end,n_strand)
#%%
from ma_mapper import gzmaf
maf_dir = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/zoonomia_241_species'
maf_file_prefix = '241-mammalian-2020v2b.maf'
maf_filepath = f'{maf_dir}/{maf_file_prefix}.{chrom}.gz'
mafindex_filepath = f'{maf_dir}/{maf_file_prefix}.{chrom}.mafindex'
target_species = 'Homo_sapiens'
maf_id = f'{target_species}.{chrom}'
index_maf = gzmaf.gzMafIndex(mafindex_filepath, maf_filepath, maf_id)
n_strand = -1 if strand == '-' else 1
fetched =index_maf.get_spliced(start,end,n_strand)
#%%
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
#%%
target_seqname = f"{target_species}.{chrom}"
expected_letters = sum(end - start for start, end in zip(start, end))
print(expected_letters)
if len(fetched) == 0:
    df = pd.DataFrame({'seqid': [target_seqname], 'seq': [Seq("N" * expected_letters)]})
if len(fetched) == 1:
    df = pd.DataFrame({'seqid': [target_seqname], 'seq': [Seq("N" * expected_letters)]})
all_seqnames = {sequence.id for multiseq in fetched for sequence in multiseq}
split_by_position = {seq_name: {} for seq_name in all_seqnames}
split_by_position
total_rec_length = 0
ref_first_strand = None
for multiseq in fetched:
    for seqrec in multiseq:
        if seqrec.id == target_seqname:
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
                    % target_seqname
                ) from None

            rec_length = len(seqrec)
            print(rec_length)
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
                "Did not find %s in alignment bundle" % (target_seqname,)
            )
    real_pos = rec_start
    edit_id = []
    edit_pos = []
    for gapped_pos in range(rec_length):
        previous_id = ''
        for seqrec in multiseq:
            
            if seqrec.id == target_seqname:
                track_val = seqrec.seq[gapped_pos]
            
            
            split_by_position[seqrec.id][real_pos] += seqrec.seq[gapped_pos]
            if previous_id == seqrec.id:
                    edit_id.append(seqrec.id)
                    edit_pos.append(real_pos)
            previous_id = seqrec.id
        if track_val != "-" and real_pos < rec_end:
            real_pos += 1
    # Debugging: Print lengths of sequences in split_by_position
    for i in range(len(edit_id)):
        _sequence=split_by_position[edit_id[i]][edit_pos[i]]
        new_sequence=process_sequence_localized(_sequence)
        split_by_position[edit_id[i]][edit_pos[i]] = new_sequence
    
    if len(split_by_position[target_seqname]) != total_rec_length:
        raise ValueError(
            "Target seqname (%s) has %s records, expected %s"
            % (
                target_seqname,
                len(split_by_position[target_seqname]),
                total_rec_length,
            )
        )

realpos_to_len = {
    pos: len(gapped_fragment)
    for pos, gapped_fragment in split_by_position[target_seqname].items()
    if len(gapped_fragment) > 1
}

seqid_list = []
seq_list = []

for seqid in all_seqnames:
    seq_split = split_by_position[seqid]
    seq_splice = []
    filler_char = "N" if seqid == target_seqname else "-"
    append = seq_splice.append
    for exonstart, exonend in zip(start, end):
        for real_pos in range(exonstart, exonend):
            if real_pos in seq_split:
                append(seq_split[real_pos])
            elif real_pos in realpos_to_len:
                append(filler_char * realpos_to_len[real_pos])
            else:
                append(filler_char)

    seqid_list.append(seqid)
    seq_list.append(Seq("".join(seq_splice))) 

if len(seq_list[seqid_list.index(target_seqname)].replace("-", "")) != expected_letters:
    raise ValueError(
        "Returning %s letters for target seqname (%s), expected %s"
        % (
            len(seq_list[seqid_list.index(target_seqname)].replace("-", "")),
            target_seqname,
            expected_letters,
        )
    )

ref_subseq_len = len(seq_list[seqid_list.index(target_seqname)])
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
#%%
def remove_hg38_gap_columns(df):
    # Find hg38 sequence
    hg38_row = df[df["seqid"].str.startswith("hg38")]
    if hg38_row.empty:
        return df  # no hg38, return unchanged
    
    hg38_seq = hg38_row.iloc[0]["seq"]
    keep_indices = [i for i, base in enumerate(hg38_seq) if base != '-']
    
    def filter_seq(seq):
        return ''.join(seq[i] for i in keep_indices)

    df = df.copy()
    df["seq"] = df["seq"].apply(filter_seq)
    return df

df_edit = remove_hg38_gap_columns(df)

# %%
