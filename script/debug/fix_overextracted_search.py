#%%
from Bio.AlignIO import MafIO
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ma_mapper import gzmaf
#%%
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
    from typing import List, Tuple

    def trim_te_fragments(positions: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
        """
        Trim overlapping transposable element (TE) fragments into non-redundant regions.

        This function is necessary for **multi-species sequence alignments** where 
        overlapping TE fragments (e.g. from RepeatMasker) can introduce duplicated sequence 
        content if not properly merged. Without this step, fragments with physical overlap 
        in one species can be duplicated in the alignment, causing:

            - False inflation of sequence length
            - Artificially high TE conservation across species
            - Misleading phylogenetic signals
            - Incorrect ancestral reconstruction

        In contrast, if working in **single-genome**, retaining overlaps 
        may preserve internal deletions or tandem repeat events and should be handled differently.

        Parameters:
        -----------
        positions : List[Tuple[int, int]]
            List of genomic (start, end) coordinates for TE fragments.

        Returns:
        --------
        List[Tuple[int, int]]
            List of non-overlapping (start, end) regions merged to cover all input spans.
        """
        if not positions:
            return []

        # Sort by start, then end
        positions.sort()
        merged = [positions[0]]

        for current_start, current_end in positions[1:]:
            last_start, last_end = merged[-1]
            if current_start <= last_end:
                # Merge overlap: take max end
                merged[-1] = (last_start, max(last_end, current_end))
            else:
                # No overlap, keep as new entry
                merged.append((current_start, current_end))

        return merged

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
    #return fetched
    expected_letters = sum(end - start for start, end in zip(starts, ends))
    #exit early if no alignment
    #print(expected_letters)
    print(fetched)
    if len(fetched) == 0:
        return pd.DataFrame({'seqid': [self._target_seqname], 'seq': [Seq("N" * expected_letters)]})

    all_seqnames = {sequence.id for multiseq in fetched for sequence in multiseq}    
    split_by_position = {seq_name: {} for seq_name in all_seqnames}
    split_by_position
    total_rec_length = 0
    ref_first_strand = None
    for multiseq in fetched:
        
        if len(multiseq) == 1 and multiseq[0].id == self._target_seqname:
            # Fast-path for human-only alignment block
            rec = multiseq[0]
            rec_start = rec.annotations["start"]
            real_pos = rec_start
            gapped_seq = rec.seq
            ref_first_strand = seqrec.annotations["strand"]
            if "-" not in gapped_seq:
                # Fully ungapped, directly assign bases by offset
                for exonstart, exonend in zip(starts, ends):
                    start_idx = exonstart - rec_start
                    end_idx = exonend - rec_start
                    for offset, base in zip(range(exonstart, exonend), gapped_seq[start_idx:end_idx]):
                        split_by_position[rec.id][offset] = base
            else:
                # Gapped version: build position map efficiently
                valid_indexes = []
                positions = []
                realpos = rec_start
                for i, base in enumerate(gapped_seq):
                    if base != "-":
                        valid_indexes.append((i, realpos))
                        positions.append(realpos)
                        realpos += 1

                for exonstart, exonend in zip(starts, ends):
                    i_start = bisect_left(positions, exonstart)
                    i_end = bisect_right(positions, exonend - 1)
                    for i in range(i_start, i_end):
                        aln_idx, real_pos = valid_indexes[i]
                        split_by_position[rec.id][real_pos] = gapped_seq[aln_idx]

            continue  # skip rest of loop for this block
        # If not fast-path, continue with your normal multi-species logic...
        for seqrec in multiseq:
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
        for i in range(len(edit_id)):
            _sequence=split_by_position[edit_id[i]][edit_pos[i]]
            new_sequence=process_sequence_localized(_sequence)
            split_by_position[edit_id[i]][edit_pos[i]] = new_sequence
        
        if len(split_by_position[self._target_seqname]) != total_rec_length:
            raise ValueError(
                "Target seqname (%s) has %s records, expected %s"
                % (
                    self._target_seqname,
                    len(split_by_position[self._target_seqname]),
                    total_rec_length,
                )
            )

    realpos_to_len = {
        pos: len(gapped_fragment)
        for pos, gapped_fragment in split_by_position[self._target_seqname].items()
        if len(gapped_fragment) > 1
    }

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
    target_index = seqid_list.index(self._target_seqname)
    if len(seq_list[target_index].replace("-", "")) != expected_letters:
        raise ValueError(
            "Returning %s letters for target seqname (%s), expected %s"
            % (
                len(seq_list[target_index].replace("-", "")),
                self._target_seqname,
                expected_letters,
            )
        )

    ref_subseq_len = len(seq_list[target_index])
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
    return df

MafIO.MafIndex.get_spliced = get_spliced_mod

#%%
def search_mod(self, starts, ends):
    # verify the provided exon coordinates
    if len(starts) != len(ends):
        raise ValueError("Every position in starts must have a match in ends")

    # Could it be safer to sort the (exonstart, exonend) pairs?
    for exonstart, exonend in zip(starts, ends):
        exonlen = exonend - exonstart
        if exonlen < 1:
            raise ValueError(
                "Exon coordinates (%d, %d) invalid: exon length (%d) < 1"
                % (exonstart, exonend, exonlen)
            )
    con = self._con

    yielded_rec_coords = set()
    # search for every exon
    for exonstart, exonend in zip(starts, ends):
        try:
            possible_bins = ", ".join(
                map(str, self._region2bin(exonstart, exonend))
            )
        except TypeError:
            raise TypeError(
                "Exon coordinates must be integers "
                "(start=%d, end=%d)" % (exonstart, exonend)
            ) from None
        result = con.execute(
            "SELECT DISTINCT start, end, offset FROM offset_data "
            "WHERE bin IN (%s) "
            "AND (end BETWEEN %s AND %s OR %s BETWEEN start AND end) "
            "ORDER BY start, end, offset ASC;"
            % (possible_bins, exonstart, exonend - 1, exonend - 1)
        )
        rows = result.fetchall()
        for rec_start, rec_end, offset in rows:
            if (rec_start, rec_end) in yielded_rec_coords:
                continue
            else:
                yielded_rec_coords.add((rec_start, rec_end))


            fetched = self._get_record(int(offset))

            for record in fetched:
                if record.id == self._target_seqname:
    
                    start = record.annotations["start"]
                    end = start + record.annotations["size"] - 1

                    if not (start == rec_start and end == rec_end):
                        raise ValueError(
                            "Expected %s-%s @ offset %s, found %s-%s"
                            % (rec_start, rec_end, offset, start, end)
                        )

            yield fetched
#%%
def search_debug(self, starts, ends):
    """Search index database for MAF records overlapping ranges provided."""
    if len(starts) != len(ends):
        raise ValueError("Every position in starts must have a match in ends")
    
    for exonstart, exonend in zip(starts, ends):
        exonlen = exonend - exonstart
        if exonlen < 1:
            raise ValueError(f"Exon coordinates ({exonstart}, {exonend}) invalid: exon length ({exonlen}) < 1")
        
        print(f"Processing exon from {exonstart} to {exonend}")
        con = self._con
        
        yielded_rec_coords = set()
        
        try:
            possible_bins = ", ".join(map(str, self._region2bin(exonstart, exonend)))
            print(f"Possible bins for start, end: {possible_bins}")
        except TypeError:
            raise TypeError(f"Exon coordinates must be integers (start={exonstart}, end={exonend})") from None
        
        sql_query = f"SELECT DISTINCT start, end, offset FROM offset_data WHERE bin IN ({possible_bins}) AND (end BETWEEN {exonstart} AND {exonend - 1} OR {exonend - 1} BETWEEN start AND end) ORDER BY start, end, offset ASC;"
        print(f"SQL query: {sql_query}")
        
        result = con.execute(sql_query)
        rows = result.fetchall()
        print(f"Rows retrieved: {rows}")
        
        for rec_start, rec_end, offset in rows:
            print(f"Processing record from {rec_start} to {rec_end} with offset {offset}")
            if (rec_start, rec_end) in yielded_rec_coords:
                continue
            else:
                yielded_rec_coords.add((rec_start, rec_end))
            
            fetched = self._get_record(int(offset))
            
            for record in fetched:
                if record.id == self._target_seqname:
                    start = record.annotations["start"]
                    end = start + record.annotations["size"] - 1
                    
                    print(f"Fetched record start: {start}, end: {end}")
                    
                    if not (start == rec_start and end == rec_end):
                        raise ValueError(f"Expected {rec_start}-{rec_end} @ offset {offset}, found {start}-{end}")
            
            print(f"Yielded record coordinates: {yielded_rec_coords}")
            yield fetched


MafIO.MafIndex.search = search_debug
#%%
def search_debug3(self, starts, ends):
    # verify the provided exon coordinates
    if len(starts) != len(ends):
        raise ValueError("Every position in starts must have a match in ends")

    # Could it be safer to sort the (exonstart, exonend) pairs?
    for exonstart, exonend in zip(starts, ends):
        exonlen = exonend - exonstart
        if exonlen < 1:
            raise ValueError(
                "Exon coordinates (%d, %d) invalid: exon length (%d) < 1"
                % (exonstart, exonend, exonlen)
            )
    con = self._con

    yielded_rec_coords = set()
    # search for every exon
    for exonstart, exonend in zip(starts, ends):
        try:
            possible_bins = ", ".join(
                map(str, self._region2bin(exonstart, exonend))
            )
        except TypeError:
            raise TypeError(
                "Exon coordinates must be integers "
                "(start=%d, end=%d)" % (exonstart, exonend)
            ) from None
        result = con.execute(
            "SELECT DISTINCT start, end, offset FROM offset_data "
            "WHERE bin IN (%s) "
            "AND (end BETWEEN %s AND %s OR %s BETWEEN start AND end) "
            "ORDER BY start, end, offset ASC;"
            % (possible_bins, exonstart, exonend - 1, exonend - 1)
        )
        rows = result.fetchall()
        for rec_start, rec_end, offset in rows:
            print(f"Processing record from {rec_start} to {rec_end} with offset {offset}")
            if (rec_start, rec_end) in yielded_rec_coords:
                continue
            else:
                yielded_rec_coords.add((rec_start, rec_end))


            fetched = self._get_record(int(offset))

            for record in fetched:
                if record.id == self._target_seqname:
    
                    start = record.annotations["start"]
                    end = start + record.annotations["size"] - 1

                    if not (start == rec_start and end == rec_end):
                        raise ValueError(
                            "Expected %s-%s @ offset %s, found %s-%s"
                            % (rec_start, rec_end, offset, start, end)
                        )

            yield fetched
MafIO.MafIndex.search = search_debug3
#%%
strand = 1
#chrom = 'chr1'
#start_list = [119563] #[156425872]
#end_list = [119944] #[156425923]
#start_list = [8000]
#end_list = [9000]
start_list = [26657850]
end_list = [26657950]
chrom = 'chrY'
maf_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/{chrom}.maf'
mafindex_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/{chrom}.mafindex'
target_species = 'hg38'
target_chrom = f'{target_species}.{chrom}'
index_maf = MafIO.MafIndex(mafindex_filepath, maf_filepath, target_chrom)
#%%
spliced_maf_full = index_maf.get_spliced(start_list, end_list, strand)
spliced_maf_full
#%%
import sqlite3

# Connect to the .mafindex file
conn = sqlite3.connect(mafindex_filepath)  # Replace with your actual path

# Create a cursor to execute queries
cursor = conn.cursor()
#%%
start_value = 26657850
end_value = 26657910
cursor.execute("SELECT * FROM offset_data WHERE start >= ? AND end <= ?;", (start_value, end_value))
rows = cursor.fetchall()
for row in rows:
    print(row)
#%%
start_list = [26657850]
end_list = [26657950]
fetched = list(index_maf.search(start_list, end_list))
print(len(fetched))  # How many alignments are returned?
for alignment in fetched:
    for seqrec in alignment:
        print(seqrec.annotations) 
    print(alignment)  # Print alignment details
#%%
expected_letters = sum(end - start for start, end in zip(start_list, end_list))
#exit early if no alignment
#print(expected_letters)
print(fetched)
if len(fetched) == 0:
    df=pd.DataFrame({'seqid': [target_chrom], 'seq': [Seq("N" * expected_letters)]})

all_seqnames = {sequence.id for multiseq in fetched for sequence in multiseq}    
split_by_position = {seq_name: {} for seq_name in all_seqnames}
split_by_position
total_rec_length = 0
ref_first_strand = None
#%%
for multiseq in fetched:
    if len(multiseq) == 1 and multiseq[0].id == target_chrom:
        # Fast-path for human-only alignment block
        rec = multiseq[0]
        rec_start = rec.annotations["start"]
        gapped_seq = rec.seq
        if ref_first_strand is None:
            print(seqrec.annotations["strand"])
            ref_first_strand = seqrec.annotations["strand"]
        if "-" not in gapped_seq:
            print('debug-homanonly: superfast route')
            # Fully ungapped, directly assign bases by offset
            for exonstart, exonend in zip(start_list, end_list):
                # Clamp the range to only what's within the alignment block
                effective_start = max(exonstart, rec_start)
                effective_end = min(exonend, rec_start + len(gapped_seq))

                if effective_start >= effective_end:
                    # Exon lies completely outside this block — skip
                    continue

                start_idx = effective_start - rec_start
                end_idx = effective_end - rec_start
                for offset, base in zip(range(effective_start, effective_end), gapped_seq[start_idx:end_idx]):
                    split_by_position[rec.id][offset] = base
        else:
            print('debug-humanonly: fast route')
            # Gapped version: build position map efficiently
            valid_indexes = []
            positions = []
            realpos = rec_start
            for i, base in enumerate(gapped_seq):
                if base != "-":
                    valid_indexes.append((i, realpos))
                    positions.append(realpos)
                    realpos += 1

            for exonstart, exonend in zip(start_list, end_list):
                i_start = bisect_left(positions, exonstart)
                i_end = bisect_right(positions, exonend - 1)
                for i in range(i_start, i_end):
                    aln_idx, real_pos = valid_indexes[i]
                    split_by_position[rec.id][real_pos] = gapped_seq[aln_idx]

        continue  # Skip rest of loop for this block
    
    for seqrec in multiseq:
        if seqrec.id == target_chrom:
            try:
                if ref_first_strand is None:
                    print(seqrec.annotations["strand"])
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
                    % target_chrom
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
                "Did not find %s in alignment bundle" % (target_chrom,)
            )
    real_pos = rec_start
    edit_id = []
    edit_pos = []
    for gapped_pos in range(rec_length):
        previous_id = ''
        for seqrec in multiseq:
            
            if seqrec.id == target_chrom:
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
    
    if len(split_by_position[target_chrom]) != total_rec_length:
        raise ValueError(
            "Target seqname (%s) has %s records, expected %s"
            % (
                target_chrom,
                len(split_by_position[target_chrom]),
                total_rec_length,
            )
        )
#%%
realpos_to_len = {
    pos: len(gapped_fragment)
    for pos, gapped_fragment in split_by_position[target_chrom].items()
    if len(gapped_fragment) > 1
}

seqid_list = []
seq_list = []

for seqid in all_seqnames:
    seq_split = split_by_position[seqid]
    seq_splice = []
    filler_char = "N" if seqid == target_chrom else "-"
    append = seq_splice.append

    for exonstart, exonend in zip(start_list, end_list):
        for real_pos in range(exonstart, exonend):
            if real_pos in seq_split:
                append(seq_split[real_pos])
            elif real_pos in realpos_to_len:
                append(filler_char * realpos_to_len[real_pos])
            else:
                append(filler_char)
    
    seqid_list.append(seqid)
    seq_list.append(Seq("".join(seq_splice))) 
target_index = seqid_list.index(target_chrom)
if len(seq_list[target_index].replace("-", "")) != expected_letters:
    raise ValueError(
        "Returning %s letters for target seqname (%s), expected %s"
        % (
            len(seq_list[target_index].replace("-", "")),
            target_chrom,
            expected_letters,
        )
    )

ref_subseq_len = len(seq_list[target_index])
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
expected_letters = sum(end - start for start, end in zip(start_list, end_list))
#exit early if no alignment
#print(expected_letters)
print(fetched)
if len(fetched) == 0:
    df=pd.DataFrame({'seqid': [target_chrom], 'seq': [Seq("N" * expected_letters)]})

all_seqnames = {sequence.id for multiseq in fetched for sequence in multiseq}    
split_by_position = {seq_name: {} for seq_name in all_seqnames}
split_by_position
total_rec_length = 0
ref_first_strand = None
# %%
from collections import Counter
from bisect import bisect_left, bisect_right
for multiseq in fetched:
    print(multiseq)
    if len(multiseq) == 1 and multiseq[0].id == target_chrom:
        # Fast-path for human-only alignment block
        rec = multiseq[0]
        rec_start = rec.annotations["start"]
        gapped_seq = rec.seq
        if ref_first_strand is None:
            ref_first_strand = seqrec.annotations["strand"]
        if "-" not in gapped_seq:
            print('debug-humanonly: superfast route')
            # Fully ungapped, directly assign bases by offset
            
            for exonstart, exonend in zip(start_list, end_list):
                # Clamp the range to only what's within the alignment block
                effective_start = max(exonstart, rec_start)
                effective_end = min(exonend, rec_start + len(gapped_seq))

                if effective_start >= effective_end:
                    # Exon lies completely outside this block — skip
                    continue

                start_idx = effective_start - rec_start
                end_idx = effective_end - rec_start
                for offset, base in zip(range(effective_start, effective_end), gapped_seq[start_idx:end_idx]):
                    split_by_position[rec.id][offset] = base
        else:
            print('debug-humanonly: fast route')
            # Gapped version: build position map efficiently
            valid_indexes = []
            positions = []
            realpos = rec_start
            for i, base in enumerate(gapped_seq):
                if base != "-":
                    valid_indexes.append((i, realpos))
                    positions.append(realpos)
                    realpos += 1

            for exonstart, exonend in zip(start_list, end_list):
                i_start = bisect_left(positions, exonstart)
                i_end = bisect_right(positions, exonend - 1)
                for i in range(i_start, i_end):
                    aln_idx, real_pos = valid_indexes[i]
                    split_by_position[rec.id][real_pos] = gapped_seq[aln_idx]

        continue  # Skip rest of loop for this block

    # Find the human reference sequence
    ref = next((r for r in multiseq if r.id == target_chrom), None)
    if ref is None:
        raise ValueError(
            "Did not find %s in alignment bundle" % (target_chrom,)
        )

    try:
        strand = ref.annotations["strand"]
    except KeyError:
        raise ValueError(
            "No strand information for target seqname (%s)" % target_chrom
        ) from None

    if ref_first_strand is None:
        ref_first_strand = strand
        if ref_first_strand not in (1, -1):
            raise ValueError("Strand must be 1 or -1")
    elif ref_first_strand != strand:
        raise ValueError(
            "Encountered strand='%s' on target seqname, expected '%s'"
            % (strand, ref_first_strand)
        )

    ref_start = ref.annotations["start"]
    ref_seq = ref.seq

    # Fast path: completely ungapped human reference
    if "-" not in ref_seq:
        print("debug: superfast route")
        for exonstart, exonend in zip(start_list, end_list):
            start_idx = exonstart - ref_start
            end_idx = exonend - ref_start
            for offset, base in zip(range(exonstart, exonend), ref_seq[start_idx:end_idx]):
                split_by_position[ref.id][offset] = base

    # Gapped path: build ref_pos_map
    print("debug: fast route")
    ref_pos_map = []
    realpos = ref_start
    for i, base in enumerate(ref_seq):
        if base != "-":
            ref_pos_map.append((i, realpos))
            realpos += 1

    ref_positions = [pos for _, pos in ref_pos_map]

    # Now extract aligned fragments for all sequences using the same alignment coordinates
    for exonstart, exonend in zip(start_list, end_list):
        i_start = bisect_left(ref_positions, exonstart)
        i_end = bisect_right(ref_positions, exonend - 1)
        if i_start == i_end:
            continue  # This exon not represented in the alignment block

        aln_indexes = [ref_pos_map[i][0] for i in range(i_start, i_end)]
        ref_coords = range(exonstart, exonstart + len(aln_indexes))

        for rec in multiseq:
            rec_id = rec.id
            gapped_seq = rec.seq
            frag = "".join(gapped_seq[i] for i in aln_indexes)

            # no need to mask as we focus on accurate extraction
            """
            if set(frag.upper()) <= {"-", "N"}:
                frag = "N" * len(frag)
            """
            for pos, base in zip(ref_coords, frag):
                split_by_position[rec_id][pos] = base
        # Count how many times each seqrec.id appears (detect duplicates)
        id_counts = Counter(seqrec.id for seqrec in multiseq)
        duplicate_ids = {seq_id for seq_id, count in id_counts.items() if count > 1}

        # Process duplicates only
        for seq_id in duplicate_ids:
            for pos, seq_str in split_by_position[seq_id].items():
                # If sequence length > 1, means concatenation happened (duplication)
                if len(seq_str) > 1:
                    new_seq = process_sequence_localized(seq_str)
                    split_by_position[seq_id][pos] = new_seq

realpos_to_len = {
    pos: len(gapped_fragment)
    for pos, gapped_fragment in split_by_position[target_chrom].items()
    if len(gapped_fragment) > 1
}

seqid_list = []
seq_list = []

for seqid in all_seqnames:
    seq_split = split_by_position[seqid]
    seq_splice = []
    filler_char = "N" if seqid == target_chrom else "-"
    append = seq_splice.append

    for exonstart, exonend in zip(start_list, end_list):
        for real_pos in range(exonstart, exonend):
            if real_pos in seq_split:
                append(seq_split[real_pos])
            elif real_pos in realpos_to_len:
                append(filler_char * realpos_to_len[real_pos])
            else:
                append(filler_char)
    
    seqid_list.append(seqid)
    seq_list.append(Seq("".join(seq_splice))) 
target_index = seqid_list.index(target_chrom)
if len(seq_list[target_index].replace("-", "")) != expected_letters:
    raise ValueError(
        "Returning %s letters for target seqname (%s), expected %s"
        % (
            len(seq_list[target_index].replace("-", "")),
            target_chrom,
            expected_letters,
        )
    )

ref_subseq_len = len(seq_list[target_index])
for seqid, seq in zip(seqid_list, seq_list):
    if len(seq) != ref_subseq_len:
        raise ValueError(
            "Returning length %s for %s, expected %s"
            % (len(seq), seqid, ref_subseq_len)
        )
#%%
# Create a DataFrame
df = pd.DataFrame({
    'seqid': seqid_list,
    'seq': [seq.reverse_complement() if strand != ref_first_strand else seq for seq in seq_list]
})
# %% REAL TE tEST
from typing import List, Tuple

def trim_overlapping_fragments(starts: List[int], ends: List[int]) -> Tuple[List[int], List[int]]:
    """
    Trim overlapping TE fragments by adjusting the start positions to avoid base duplication.

    This is used to **prevent per-base duplication** when processing overlapping TE fragments,
    especially in multi-species alignments. It ensures that each base in the reference genome
    is only used once, even if it is part of multiple fragment calls.

    This function is necessary for **multi-species sequence alignments** where 
        overlapping TE fragments (e.g. from RepeatMasker) can introduce duplicated sequence 
        content if not properly merged. Without this step, fragments with physical overlap 
        in one species can be duplicated in the alignment, causing:

            - False inflation of sequence length
            - Artificially high TE conservation across species
            - Misleading phylogenetic signals
            - Incorrect ancestral reconstruction

        In contrast, if working in **single-genome**, retaining overlaps 
        may preserve internal deletions or tandem repeat events and should be handled differently.

    Parameters:
    -----------
    starts : List[int]
        List of start positions.
    ends : List[int]
        List of end positions (must be same length as `starts`).

    Returns:
    --------
    Tuple[List[int], List[int]]
        Adjusted (trimmed) start and end positions with no overlapping bases.
    """
    if not starts or not ends or len(starts) != len(ends):
        raise ValueError("Starts and ends must be non-empty and of equal length.")

    trimmed_starts = [starts[0]]
    trimmed_ends = [ends[0]]

    for i in range(1, len(starts)):
        prev_end = trimmed_ends[-1]
        curr_start = starts[i]
        curr_end = ends[i]

        # Trim current start to be >= previous end
        adjusted_start = max(curr_start, prev_end)
        if adjusted_start >= curr_end:
            # Skip completely overlapped region
            continue
        trimmed_starts.append(adjusted_start)
        trimmed_ends.append(curr_end)

    return trimmed_starts, trimmed_ends

        
strand = -1
#start_list = [81195288]
#end_list = [81195660]
start_list = [111584863, 111584915, 111585238]
end_list = [111584927, 111584958, 111585407]
chrom = 'chr13'
maf_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/{chrom}.maf'
mafindex_filepath = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/{chrom}.mafindex'
target_species = 'hg38'
target_chrom = f'{target_species}.{chrom}'
index_maf = MafIO.MafIndex(mafindex_filepath, maf_filepath, target_chrom)
#%%
start_list,end_list=trim_overlapping_fragments(start_list, end_list)
fetched = list(index_maf.search(start_list, end_list))
print(len(fetched))  # How many alignments are returned?
#for alignment in fetched:
#    for seqrec in alignment:
#        print(seqrec.annotations) 
#    print(alignment)  # Print alignment details
expected_letters = sum(end - start for start, end in zip(start_list, end_list))
#exit early if no alignment
#print(expected_letters)
print(fetched)
if len(fetched) == 0:
    df=pd.DataFrame({'seqid': [target_chrom], 'seq': [Seq("N" * expected_letters)]})

all_seqnames = {sequence.id for multiseq in fetched for sequence in multiseq}    
split_by_position = {seq_name: {} for seq_name in all_seqnames}
#split_by_position
total_rec_length = 0
ref_first_strand = None
#%%
from collections import Counter, defaultdict
from bisect import bisect_left, bisect_right
split_by_position = defaultdict(lambda: defaultdict(str))
#Precompute all exon positions (once before loop)
exon_positions = []
for start, end in zip(start_list, end_list):
    exon_positions.extend(range(start, end))
exon_positions = sorted(set(exon_positions))  # Ensure unique & sorted
for multiseq in fetched:
    print(multiseq)
    if len(multiseq) == 1 and multiseq[0].id == target_chrom:
        # Fast-path for human-only alignment block
        rec = multiseq[0]
        rec_start = rec.annotations["start"]
        gapped_seq = rec.seq
        if ref_first_strand is None:
            ref_first_strand = rec.annotations["strand"]
        if "-" not in gapped_seq:
            #print('debug-humanonly: superfast route')
            # Fully ungapped, directly assign bases by offset
            
            for exonstart, exonend in zip(start_list, end_list):
                # Clamp the range to only what's within the alignment block
                effective_start = max(exonstart, rec_start)
                effective_end = min(exonend, rec_start + len(gapped_seq))
                print(effective_start, effective_end)
                if effective_start >= effective_end:
                    # Exon lies completely outside this block — skip
                    continue

                start_idx = effective_start - rec_start
                end_idx = effective_end - rec_start
                for offset, base in zip(range(effective_start, effective_end), gapped_seq[start_idx:end_idx]):
                    split_by_position[rec.id][offset] = base
        else:
            #print('debug-humanonly: fast route')
            # Gapped version: build position map efficiently
            valid_indexes = []
            positions = []
            realpos = rec_start
            for i, base in enumerate(gapped_seq):
                if base != "-":
                    valid_indexes.append((i, realpos))
                    positions.append(realpos)
                    realpos += 1

            for exonstart, exonend in zip(start_list, end_list):
                i_start = bisect_left(positions, exonstart)
                i_end = bisect_right(positions, exonend - 1)
                for i in range(i_start, i_end):
                    aln_idx, real_pos = valid_indexes[i]
                    print(aln_idx, real_pos)
                    split_by_position[rec.id][real_pos] += gapped_seq[aln_idx]

        continue  # Skip rest of loop for this block

    # Find the human reference sequence
    ref = next((r for r in multiseq if r.id == target_chrom), None)
    if ref is None:
        raise ValueError(
            "Did not find %s in alignment bundle" % (target_chrom,)
        )
    
    try:
        if ref_first_strand is None:
            ref_first_strand = ref.annotations["strand"]
            if ref_first_strand not in (1, -1):
                raise ValueError("Strand must be 1 or -1")
        elif ref_first_strand != ref.annotations["strand"]:
            raise ValueError(
                "Encountered strand='%s' on target seqname, expected '%s'"
                % (strand, ref_first_strand)
            )
    except KeyError:
        raise ValueError(
            "No strand information for target seqname (%s)" % target_chrom
        ) from None

    ref_start = ref.annotations["start"]
    ref_seq = ref.seq

    # Gapped path: build ref_pos_map
    #print("debug-humanref: fast route")
    ref_pos_map = []
    realpos = ref_start
    for i, base in enumerate(ref_seq):
        if base != "-":
            ref_pos_map.append((i, realpos))
            realpos += 1
    
    block_positions  = [pos for _, pos in ref_pos_map]
    """
    for i in ref_positions:
        print(i)
    """
    # Now extract aligned fragments for all sequences using the same alignment coordinates
    for exonstart, exonend in zip(start_list, end_list):
        
        i_start = bisect_left(block_positions, exonstart)
        i_end = bisect_right(block_positions, exonend - 1)
        print(exonstart,exonend,'after bisect',i_start,i_end)
        if i_start == i_end:
            continue  # This exon not represented in the alignment block

        aln_indexes = [ref_pos_map[i][0] for i in range(i_start, i_end)]
        ref_coords = block_positions[i_start:i_end]  # genomic coords covered in this exon by this block
        """
        print(len(block_positions))
        print(block_positions)
        print(aln_indexes)
        print(ref_coords)
        """
        for rec in multiseq:
            rec_id = rec.id
            gapped_seq = rec.seq
            frag = "".join(gapped_seq[i] for i in aln_indexes)

            # no need to mask as we focus on accurate extraction
            """
            if set(frag.upper()) <= {"-", "N"}:
                frag = "N" * len(frag)
            """
            #print(ref_coords, frag)
            for pos, base in zip(ref_coords, frag):
                
                split_by_position[rec_id][pos] += base
                if rec_id == target_chrom:
                    print(pos,split_by_position[rec_id][pos])
        # Count how many times each seqrec.id appears (detect duplicates)
        id_counts = Counter(seqrec.id for seqrec in multiseq)
        duplicate_ids = {seq_id for seq_id, count in id_counts.items() if count > 1}

        # Process duplicates only
        for seq_id in duplicate_ids:
            for pos, seq_str in split_by_position[seq_id].items():
                # If sequence length > 1, means concatenation happened (duplication)
                if len(seq_str) > 1:
                    new_seq = process_sequence_localized(seq_str)
                    split_by_position[seq_id][pos] = new_seq
realpos_to_len = {
    pos: len(gapped_fragment)
    for pos, gapped_fragment in split_by_position[target_chrom].items()
    if len(gapped_fragment) > 1
}

seqid_list = []
seq_list = []

for seqid in all_seqnames:
    seq_split = split_by_position[seqid]
    seq_splice = []
    filler_char = "N" if seqid == target_chrom else "-"
    append = seq_splice.append

    for exonstart, exonend in zip(start_list, end_list):
        for real_pos in range(exonstart, exonend):
            if real_pos in seq_split:
                append(seq_split[real_pos])
            elif real_pos in realpos_to_len:
                append(filler_char * realpos_to_len[real_pos])
            else:
                append(filler_char)
    
    seqid_list.append(seqid)
    seq_list.append(Seq("".join(seq_splice))) 
target_index = seqid_list.index(target_chrom)
if len(seq_list[target_index].replace("-", "")) != expected_letters:
    raise ValueError(
        "Returning %s letters for target seqname (%s), expected %s"
        % (
            len(seq_list[target_index].replace("-", "")),
            target_chrom,
            expected_letters,
        )
    )

ref_subseq_len = len(seq_list[target_index])
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
df
# %%
