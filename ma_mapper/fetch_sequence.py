#%%
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from concurrent.futures import ProcessPoolExecutor
#%%
######################run directly from windows
#def parallel_init(metadata_var, records_var):
#    global metadata
#    metadata = metadata_var
#    global records 
#    records = records_var

def extract_sequence(genoname, genostart, genoend, strand):
    #seqname = '::'.join([meta_id,genoname,str(genostart),str(genoend),strand])
    chromosome_extract=records[genoname]
    if strand == '+':
        seq_string = str(chromosome_extract[genostart:genoend].seq)
    else:
        seq_string = str(chromosome_extract[genostart:genoend].seq.reverse_complement())
    return seq_string
    #seq_record = SeqRecord(Seq(''.join(seq_string)),seqname , '', '')
    #return seq_record

def fetch_sequence(metadata_input,source_fasta,output_filepath = None, save_to_file=False, custom_id = False):
    global metadata
    import os
    if isinstance(metadata_input, str):
        if (os.path.isfile(metadata_input) == True):
            metadata = pd.read_csv(metadata_input, delim_whitespace=True)
        else:
            print('metadata file not found')
    else:
        metadata = metadata_input

    if output_filepath is None:
        if isinstance(metadata_input, str):
            output_dir = '/'.join(str.split(metadata_input, sep ='/')[:-1]) 
        else:
            output_dir = os.path.dirname(os.path.abspath(__file__))
    else:
        output_dir = '/'.join(str.split(output_filepath, sep ='/')[:-1]) 
    import logging
    log_path = output_dir+'/'+'seq_extract.log'
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
    seq_records = []
    
    if custom_id == False:
        meta_id = 'sample_n'+ metadata.index.astype(str)
        metadata['meta_id'] = meta_id
    else:

        metadata['meta_id'] = metadata.iloc[:,4]
        meta_id = metadata.iloc[:,4].unique()
    ############### run directly from windows    with ProcessPoolExecutor(max_workers=nthread, initializer=parallel_init, initargs= (metadata, records)) as executor:
    #with ProcessPoolExecutor(max_workers=nthread) as executor:
    #    results = list(executor.map(extract_sequence, chrom, start, end, strand, meta_id))
    #for result in results:
    #    seq_records.append(result)
    for uniq_meta_id in meta_id:
        metadata_by_id = metadata[metadata.meta_id == uniq_meta_id]
        seq_strings = []
        for idx, row in metadata_by_id.iterrows():       
            chrom = row.iloc[0]
            start = row.iloc[1]
            end = row.iloc[2]
            strand = row.strand
            seq_string = extract_sequence(chrom, start, end, strand)
            seq_strings.append(seq_string)
        seqname = '::'.join([uniq_meta_id,chrom,str(min(metadata_by_id.iloc[:,1])),str(max(metadata_by_id.iloc[:,2])),strand])
        seq_record = SeqRecord(Seq(''.join(seq_strings)),seqname , '', '')
        seq_records.append(seq_record)
    if save_to_file == True:
        with open(output_filepath, "w") as output_handle:
            SeqIO.write(seq_records, output_handle, "fasta")
            logging.info('done, saving te sequences at: '+output_filepath)
    else:
        logging.info('done, returning seq_records as object')
        return seq_records
#%%
def main():
    source_fasta = '/home/pc575/phd_project_development/data/hg38_fasta/hg38.fa'
    metadata = '/home/pc575/phd_project_development/data/_ma_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/mer11a_coord_with_id.txt'
    output = '/home/pc575/phd_project_development/data/_ma_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/custom_id.fasta'
    fetch_sequence(metadata,source_fasta, output_filepath = output, custom_id= True)
#%%
if __name__ == '__main__':
    main()
#%%
