#%%
from tarfile import RECORDSIZE
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from concurrent.futures import ProcessPoolExecutor
import logging
import os
#global metadata
#global records 
#records = {}
#metadata = None
#%%
######################run directly from windows
#def parallel_init(metadata_var, records_var):
#    global metadata
#    metadata = metadata_var
#    global records 
#    records = records_var

def extract_sequence(metadata_row):
    row=metadata.iloc[metadata_row]
    genoname = row.genoName
    genostart = row.genoStart
    genoend = row.genoEnd
    strand = row.strand
    seqname = '::'.join([row.id,genoname,str(genostart),str(genoend),strand])
    chromosome_extract=records[genoname]
    if strand == '+':
        seq_string = str(chromosome_extract[genostart:genoend].seq)
    else:
        seq_string = str(chromosome_extract[genostart:genoend].seq.reverse_complement())
    seq_record = SeqRecord(Seq(''.join(seq_string)),seqname , '', '')
    return seq_record

def fetch_sequence(metadata_input,source_fasta,output_filepath = None, nthread = 1):
    global metadata
    if (os.path.isfile(metadata_input) == True):
        metadata = pd.read_csv(metadata_input, delim_whitespace=True)
    else:
        metadata = metadata_input
    if output_filepath is None:
        output_filepath = '/'.join(str.split(metadata_input, sep ='/')[:-1]) + '/seqrecords.fasta'
    log_path = output_filepath+'.log'
    #setup logger
    logging.root.handlers = []
    logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    handlers = [
                        logging.FileHandler(log_path, mode = 'a'),
                        logging.StreamHandler()
                        ]
                    )
    logging.info('read sequence source: '+ source_fasta)
    global records
    records = SeqIO.to_dict(SeqIO.parse(open(source_fasta), 'fasta'))
    seq_records = list()
    metadata_rows = metadata.index.to_list()
############### run directly from windows    with ProcessPoolExecutor(max_workers=nthread, initializer=parallel_init, initargs= (metadata, records)) as executor:
    with ProcessPoolExecutor(max_workers=nthread) as executor:
        results = list(executor.map(extract_sequence, metadata_rows))
    for result in results:
        seq_records.append(result)
    with open(output_filepath, "w") as output_handle:
        SeqIO.write(seq_records, output_handle, "fasta")
        logging.info('done, saving te sequences at: '+output_filepath)   
#%%
def main():
    source_fasta = '/home/pc575/phd_project_development/data/hg38_fasta/hg38.fa'
    metadata = '/home/pc575/phd_project_development/data/ma_mapper_output/mer11a_coord.txt'
    fetch_sequence(metadata,source_fasta, nthread= 2)
#%%
if __name__ == '__main__':
    main()
#%%