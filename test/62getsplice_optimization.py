import pandas as pd
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
    def process_row(row):
        sequence = row['sequence'].upper()
        if len(sequence) > 1:
            base_counts = Counter(sequence)
            max_count = max(base_counts.values())
            consensus_bases = [base for base, count in base_counts.items() if count == max_count]
            return convert_to_iupac(consensus_bases)
        else:
            return sequence.upper()

    data['consensus'] = data.apply(process_row, axis=1)
    return data

def get_spliced_mod(self, starts, ends, strand=1):
    if strand not in (1, -1):
        raise ValueError("Strand must be 1 or -1, got %s" % strand)

    fetched = list(self.search(starts, ends))
    expected_letters = sum(end - start for start, end in zip(starts, ends))

    if len(fetched) == 0:
        return [SeqRecord(Seq("N" * expected_letters), id=self._target_seqname)]

    all_seqnames = {sequence.id for multiseq in fetched for sequence in multiseq}
    split_by_position = pd.DataFrame(columns=['seqname', 'position', 'sequence'])

    total_rec_length = 0
    ref_first_strand = None

    for multiseq in fetched:
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
                            "expected '%s'" % (seqrec.annotations["strand"], ref_first_strand)
                        )
                except KeyError:
                    raise ValueError(
                        "No strand information for target seqname (%s)" % self._target_seqname
                    ) from None

                rec_length = len(seqrec)
                rec_start = seqrec.annotations["start"]
                ungapped_length = seqrec.annotations["size"]
                rec_end = rec_start + ungapped_length - 1
                total_rec_length += ungapped_length

                for seqrec in multiseq:
                    new_rows = pd.DataFrame({
                        'seqname': seqrec.id,
                        'position': list(range(rec_start, rec_end + 1)),
                        'sequence': [""] * (rec_end - rec_start + 1)
                    })
                    split_by_position = pd.concat([split_by_position, new_rows], ignore_index=True)

                break
            else:
                raise ValueError(
                    "Did not find %s in alignment bundle" % (self._target_seqname,)
                )

        real_pos = rec_start
        for gapped_pos in range(rec_length):
            for seqrec in multiseq:
                if seqrec.id == self._target_seqname:
                    track_val = seqrec.seq[gapped_pos]
                split_by_position.loc[
                    (split_by_position['seqname'] == seqrec.id) & (split_by_position['position'] == real_pos),
                    'sequence'
                ] += seqrec.seq[gapped_pos]
            if track_val != "-" and real_pos < rec_end:
                real_pos += 1

    # Process the sequences with the new function
    split_by_position = process_sequences_work(split_by_position)

    if len(split_by_position[split_by_position['seqname'] == self._target_seqname]) != total_rec_length:
        raise ValueError(
            "Target seqname (%s) has %s records, expected %s" % (
                self._target_seqname,
                len(split_by_position[split_by_position['seqname'] == self._target_seqname]),
                total_rec_length,
            )
        )

    # Splice together the exons
    subseq = {}
    for seqid in all_seqnames:
        seq_split = split_by_position[split_by_position['seqname'] == seqid]
        seq_splice = []
        filler_char = "N" if seqid == self._target_seqname else "-"
        for exonstart, exonend in zip(starts, ends):
            for real_pos in range(exonstart, exonend):
                seq_row = seq_split[seq_split['position'] == real_pos]
                if not seq_row.empty:
                    seq_splice.append(seq_row['consensus'].values[0])
                else:
                    seq_splice.append(filler_char)
        subseq[seqid] = "".join(seq_splice)

    if len(subseq[self._target_seqname].replace("-", "")) != expected_letters:
        raise ValueError(
            "Returning %s letters for target seqname (%s), expected %s" % (
                len(subseq[self._target_seqname].replace("-", "")),
                self._target_seqname,
                expected_letters,
            )
        )

    ref_subseq_len = len(subseq[self._target_seqname])
    for seqid, seq in subseq.items():
        if len(seq) != ref_subseq_len:
            raise ValueError(
                "Returning length %s for %s, expected %s" % (len(seq), seqid, ref_subseq_len)
            )

    result_multiseq = []
    for seqid, seq in subseq.items():
        seq = Seq(seq)
        seq = seq if strand == ref_first_strand else seq.reverse_complement()
        result_multiseq.append(SeqRecord(seq, id=seqid, name=seqid, description=""))

    return result_multiseq