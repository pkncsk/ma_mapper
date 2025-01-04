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

    
