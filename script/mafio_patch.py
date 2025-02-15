#%%
import concurrent
from itertools import islice, repeat
from ma_mapper import extract_maf
import importlib
importlib.reload(extract_maf)
from Bio.AlignIO import MafIO
try:
    from sqlite3 import dbapi2
except ImportError:
    dbapi2 = None
MAFINDEX_VERSION = 2
import time
#%%
import io
import subprocess
from subprocess import Popen, PIPE
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

    def _MafIndex__make_new_index(self):
        """Read MAF file and generate SQLite index (PRIVATE)."""
        # make the tables
        self._con.execute("CREATE TABLE meta_data (key TEXT, value TEXT);")
        self._con.execute(
            "INSERT INTO meta_data (key, value) VALUES (?, ?);",
            ("version", MAFINDEX_VERSION),
        )
        self._con.execute(
            "INSERT INTO meta_data (key, value) VALUES ('record_count', -1);"
        )
        self._con.execute(
            "INSERT INTO meta_data (key, value) VALUES (?, ?);",
            ("target_seqname", self._target_seqname),
        )
        # Determine whether to store maf file as relative to the index or absolute
        if not os.path.isabs(self._maf_file) and not os.path.isabs(self._index_filename):
            mafpath = os.path.relpath(self._maf_file, self._relative_path).replace(
                os.path.sep, "/"
            )
        elif (
            os.path.dirname(os.path.abspath(self._maf_file)) + os.path.sep
        ).startswith(self._relative_path + os.path.sep):
            mafpath = os.path.relpath(self._maf_file, self._relative_path).replace(
                os.path.sep, "/"
            )
        else:
            mafpath = os.path.abspath(self._maf_file)
        self._con.execute(
            "INSERT INTO meta_data (key, value) VALUES (?, ?);",
            ("filename", mafpath),
        )
        
        # Add length_to_next column for precomputed lengths (initialized as NULL)
        self._con.execute(
            "CREATE TABLE offset_data (bin INTEGER, start INTEGER, end INTEGER, offset INTEGER, length_to_next INTEGER DEFAULT NULL);"
        )

        insert_count = 0

        # iterate over the entire file and insert in batches
        mafindex_func = self._MafIndex__maf_indexer()

        while True:
            batch = list(islice(mafindex_func, 100))
            if not batch:
                break

            # Insert batch data without length_to_next for now (defaults to NULL)
            self._con.executemany(
                "INSERT INTO offset_data (bin, start, end, offset) VALUES (?,?,?,?);",
                batch,
            )
            self._con.commit()
            insert_count += len(batch)

        # Update length_to_next column for precomputed lengths
        self._con.execute(
            """
            WITH cte AS (
                SELECT offset, 
                    LEAD(offset) OVER (ORDER BY offset ASC) - offset AS length_to_next
                FROM offset_data
            )
            UPDATE offset_data
            SET length_to_next = cte.length_to_next
            FROM cte
            WHERE offset_data.offset = cte.offset;
            """
        )

        # then make indexes on the relevant fields
        self._con.execute("CREATE INDEX IF NOT EXISTS bin_index ON offset_data(bin);")
        self._con.execute(
            "CREATE INDEX IF NOT EXISTS start_index ON offset_data(start);"
        )
        self._con.execute("CREATE INDEX IF NOT EXISTS end_index ON offset_data(end);")

        self._con.execute(
            f"UPDATE meta_data SET value = '{insert_count}' WHERE key = 'record_count'"
        )

        self._con.commit()

        return insert_count
    @staticmethod
    def seek_and_read(offset, length_to_next, maf_file):
        
        #start_seek = time.time()    
        """
        Seek to a specific offset in the compressed file and read the block.
        Query the next offset from SQLite
        start_sql = time.time()    
        cursor = self._con.cursor()
        cursor.execute(
            "SELECT offset FROM offset_data WHERE offset > ? ORDER BY offset ASC LIMIT 1;",
            (offset,),
        )
        next_offset_row = cursor.fetchone()

        Determine the length of the block
        if next_offset_row:
            next_offset = next_offset_row[0]
            print(f'next offset: {next_offset}')
            length = next_offset - offset
        else:
        # If no next offset, read until EOF
            length = None  # Indicates to read to the end or a large chunk
        print("sql time:", time.time() - start_sql)
        """
        #start_gz = time.time()    
        data = b""
        if length_to_next is not None:
            cmd = [
                gztool_path,
                "-b", str(offset),  # Offset to seek
                "-r", str(length_to_next),  # Byte length to read
                maf_file,
            ]
        else:
            cmd = [
                gztool_path,
                "-b", str(offset),  # Offset to seek
                maf_file,
            ]
        #print(f'gztool time {time.time()-start_gz}')
        
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        #start_read = time.time()
        data += result.stdout
        #print(f'read time {time.time()-start_read}')
        start_io = time.time()
        from io import StringIO
        handle = StringIO(data.decode())
        #print(f'stringio time {time.time()-start_io}')
        try:
            #start_parse = time.time() 
            record = next(MafIO.MafIterator(handle))
            #print("parse time:", time.time() - start_parse)
            #print("seek time:", time.time() - start_seek)
            return record
        except ValueError:
            raise RuntimeError("Failed to parse block.")
        """
        
        with Popen(cmd, stdout=PIPE, stderr=PIPE, text = True) as proc:
            try:
                #start_parse = time.time() 
                record = next(MafIO.MafIterator(proc.stdout))
                #print("parse time:", time.time() - start_parse)
                #print("seek time:", time.time() - start_seek)
                return record
            except ValueError:
                raise RuntimeError("Failed to parse block.")
        
        
        with Popen(cmd, stdout=PIPE, stderr=PIPE) as proc:
        # Capture gztool's output into a BytesIO stream
            
            #start_read = time.time()
            raw_data = proc.stdout.read()
            #print(f'read time {time.time()-start_read}')
            
            #start_io = time.time()
            handle = io.StringIO(raw_data.decode())
            #print(f'stringio time {time.time()-start_io}')
            
            try:
                start_parse = time.time() 
                record = next(MafIO.MafIterator(handle))
                #print("parse time:", time.time() - start_parse)
                #print("seek time:", time.time() - start_seek)
                return record
            except ValueError:
                raise RuntimeError("Failed to parse block.")
        """
    @staticmethod
    def fetch_align(row,maf_file, target_seqname):
        rec_start, rec_end, offset, length_to_next=row
        fetched = gzMafIndex.seek_and_read(int(offset), int(length_to_next), maf_file)

        for record in fetched:
            if record.id == target_seqname:
                start = record.annotations["start"]
                end = start + record.annotations["size"] - 1

                if not (start == rec_start and end == rec_end):
                    raise ValueError(
                        "Expected %s-%s @ offset %s, found %s-%s"
                        % (rec_start, rec_end, offset, start, end)
                    )
        return fetched

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
            #start_sql = time.time()
            result = con.execute(
                """
                SELECT DISTINCT start, 
                   end, 
                   offset, 
                   length_to_next
                FROM offset_data
                WHERE bin IN (%s)
                AND (end BETWEEN %s AND %s OR %s BETWEEN start AND end)
                ORDER BY start, end, offset ASC;
                """ % (possible_bins, exonstart, exonend - 1, exonend - 1)
            )
            
            rows = result.fetchall()
            #print(f'sql time:{time.time()-start_sql}')
            print(f'unfiltered block {len(rows)}')
            #filter_start = time.time()
            filtered_count = 0
            filtered_rows = []
            total_length = 0
            for rec_start, rec_end, offset, length_to_next in rows:
                """
                print(f'recs info: {rec_start}\t{rec_end}')
                print(f'offsets: {offset},{length_to_next}')
                """
                if (rec_start, rec_end) in yielded_rec_coords:
                    filtered_count = filtered_count+1
                    continue
                else:
                    yielded_rec_coords.add((rec_start, rec_end))
                    filtered_rows.append((rec_start, rec_end, offset, length_to_next))
            
            #print(f'filter time {time.time()-filter_start}')
            #start_fetch = time.time()
            """
            for row in filtered_rows:
                fetched = fetch_align(row, self._maf_file, self._target_seqname)
                yield fetched
            """
            with concurrent.futures.ProcessPoolExecutor(max_workers = 60) as executor:
                results=executor.map(gzMafIndex.fetch_align, filtered_rows, repeat(self._maf_file), repeat(self._target_seqname))
            
            for result in results:
                #print(f"Alignment with {len(result)} rows and {result.get_alignment_length()} columns")
                total_length = total_length+result.get_alignment_length()
                yield result

            #print(f'fetch time {time.time()-start_fetch}')
            print(f'filtered {filtered_count}')
            print(f'total length {total_length}')


#%% test mafIO
start_total = time.time()  

chrom = 'chr1'
target_species = 'hg38'
maf_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/chr1.maf.gz'
mafindex_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/chr1.mafindex'
"""
target_species = 'Homo_sapiens'
maf_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/zoonomia_241_species/241-mammalian-2020v2b.maf.chr1.gz'
mafindex_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/zoonomia_241_species/241-mammalian-2020v2b.maf.chr1.mafindex'
"""
strand = '-'
start = 4381114
end = 4382190
maf_id = f'{target_species}.{chrom}'
    
index_maf = gzMafIndex(mafindex_file, maf_file, maf_id) 
n_strand = -1 if strand == '-' else 1
print([start],[end],n_strand)
results =index_maf.get_spliced([start],[end],n_strand)
print("total time:", time.time() - start_total)
 #%%
start_total = time.time()  
start_list = [start]
end_list = [end]
start_flanked=[min(start_list)-5000] + start_list + [max(end_list)]
end_flanked = [min(start_list)] + end_list + [max(end_list)+5000]
results =index_maf.get_spliced(start_flanked,end_flanked,n_strand)
print("total time:", time.time() - start_total)
#%% old
target_species = 'hg38'
#target_species = 'Homo_sapiens'
chrom = 'chr2'
maf_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/chr2.maf'
mafindex_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/chr2.mafindex'
#maf_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/zoonomia_241_species/241-mammalian-2020v2b.maf.chr1'
#mafindex_file = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/zoonomia_241_species/241-mammalian-2020v2b.maf.chr1.mafindex'

strand = '-'
start = 4381114
end = 4382190
maf_id = f'{target_species}.{chrom}'
    
index_maf = MafIO.MafIndex(mafindex_file, maf_file, maf_id) 
n_strand = -1 if strand == '-' else 1
results_old =index_maf.get_spliced([start],[end],n_strand)
#%%
start_list = [start]
end_list = [end]
start_flanked=[min(start_list)-5000] + start_list + [max(end_list)]
end_flanked = [min(start_list)] + end_list + [max(end_list)+5000]
results =index_maf.get_spliced(start_flanked,end_flanked,n_strand)
#%%
import subprocess

def is_iterable(obj):
    """Check if an object is iterable."""
    try:
        iter(obj)
        return True
    except TypeError:
        return False

# Run gztool and capture stdout
cmd = ["/rds/project/rds-XrHDlpCeVDg/users/pakkanan/dev/installation/gztool/gztool", "-b", "84252", "-r", "200", "/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/chrY.maf.gz"]
result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

# Get the output as bytes
stdout_data = result.stdout
# %%
import subprocess
from Bio.AlignIO.MafIO import MafIterator

cmd = ["/rds/project/rds-XrHDlpCeVDg/users/pakkanan/dev/installation/gztool/gztool", "-b", "31", "-r", "84221", "/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447/chrY.maf.gz"]

try:
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    # Decode and prepare the output for MafIterator
    output = result.stdout.decode()  # Assuming UTF-8 encoding
    print("Captured Output:\n", output)  # Print to inspect format

    # Test MafIterator with the output
    from io import StringIO
    handle = StringIO(output)
    for record in MafIterator(handle):
        print(record)  # Inspect the parsed record
except subprocess.CalledProcessError as e:
    print("Error:", e.stderr.decode())
# %%
