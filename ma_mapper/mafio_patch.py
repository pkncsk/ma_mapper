#%%
from ma_mapper import extract_maf
import importlib
import numpy as np
import pandas as pd
importlib.reload(extract_maf)
from Bio.AlignIO import MafIO
try:
    from sqlite3 import dbapi2
except ImportError:
    dbapi2 = None
MAFINDEX_VERSION = 2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#%%
import subprocess
gztool_path = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/dev/installation/gztool/gztool'
#%%
import os
class gzMafIndex(MafIO.MafIndex):
    def _MafIndex__check_existing_db(self):
        try:
            idx_version = int(
                self._con.execute(
                    "SELECT value FROM meta_data WHERE key = 'version'"
                ).fetchone()[0]
            )
            if idx_version != MAFINDEX_VERSION:
                msg = "\n".join(
                    [
                        "Index version (%s) incompatible with this version "
                        "of MafIndex" % idx_version,
                        "You might erase the existing index %s "
                        "for it to be rebuilt." % self._index_filename,
                    ]
                )
                raise ValueError(msg)

            filename = self._con.execute(
                "SELECT value FROM meta_data WHERE key = 'filename'"
            ).fetchone()[0]
            if os.path.isabs(filename):
                tmp_mafpath = filename
            else:
                tmp_mafpath = os.path.join(
                    self._relative_path, filename.replace("/", os.path.sep)
                )
            #amend .gz filetype
            if os.path.splitext(self._maf_file)[1] == ".gz":
                normalized_path = os.path.splitext(self._maf_file)[0]
            else:
                normalized_path = self._maf_file
            if tmp_mafpath != os.path.abspath(normalized_path):
                raise ValueError(
                    f"Index uses a different file ({filename} != {normalized_path})"
                )

            db_target = self._con.execute(
                "SELECT value FROM meta_data WHERE key = 'target_seqname'"
            ).fetchone()[0]
            if db_target != self._target_seqname:
                raise ValueError(
                    "Provided database indexed for %s, expected %s"
                    % (db_target, self._target_seqname)
                )

            record_count = int(
                self._con.execute(
                    "SELECT value FROM meta_data WHERE key = 'record_count'"
                ).fetchone()[0]
            )
            if record_count == -1:
                raise ValueError("Unfinished/partial database provided")

            records_found = int(
                self._con.execute("SELECT COUNT(*) FROM offset_data").fetchone()[0]
            )
            if records_found != record_count:
                raise ValueError(
                    "Expected %s records, found %s.  Corrupt index?"
                    % (record_count, records_found)
                )

            return records_found

        except (dbapi2.OperationalError, dbapi2.DatabaseError) as err:
            raise ValueError(f"Problem with SQLite database: {err}") from None
        
    def seek_and_read(self, offset):
        """Seek to a specific offset in the compressed file and read the block."""
        # Query the next offset from SQLite
        cursor = self._con.cursor()
        cursor.execute(
            "SELECT offset FROM offset_data WHERE offset > ? ORDER BY offset ASC LIMIT 1;",
            (offset,),
        )
        next_offset_row = cursor.fetchone()

        # Determine the length of the block
        if next_offset_row:
            next_offset = next_offset_row[0]
            #print(next_offset)
            length = next_offset - offset
        else:
            # If no next offset, read until EOF
            length = None  # Indicates to read to the end or a large chunk

        data = b""
        if length is not None:
            cmd = [
                gztool_path,
                "-b", str(offset),  # Offset to seek
                "-r", str(length),  # Byte length to read
                self._maf_file,
            ]
        else:
            cmd = [
                gztool_path,
                "-b", str(offset),  # Offset to seek
                self._maf_file,
            ]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        data += result.stdout

        from io import StringIO
        handle = StringIO(data.decode())
        try:
            record = next(MafIO.MafIterator(handle))
            return record
        except ValueError:
            raise RuntimeError("Failed to parse block.")
        
    def search(self, starts, ends):
        if len(starts) != len(ends):
            raise ValueError("Every position in starts must have a match in ends")
        for exonstart, exonend in zip(starts, ends):
            exonlen = exonend - exonstart
            if exonlen < 1:
                raise ValueError(
                    "Exon coordinates (%d, %d) invalid: exon length (%d) < 1"
                    % (exonstart, exonend, exonlen)
                )
        con = self._con
        yielded_rec_coords = set()
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
                #print(f'{rec_start}\t{rec_end}\t{offset}')
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

    def _get_record(self, offset):
        if os.path.splitext(self._maf_file)[1] == ".gz":
            record = self.seek_and_read(offset)
        else:
            self._maf_fp.seek(offset)
            record = next(self._mafiter)
        return record

    def get_spliced(self, starts, ends, strand=1):
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
        #return fetched
        expected_letters = sum(end - start for start, end in zip(starts, ends))
        if len(fetched) == 0:
            return [SeqRecord(Seq("N" * expected_letters), id=self._target_seqname)]
        all_seqnames = {sequence.id for multiseq in fetched for sequence in multiseq}
        split_by_position = {seq_name: {} for seq_name in all_seqnames}
        split_by_position
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
        return df
