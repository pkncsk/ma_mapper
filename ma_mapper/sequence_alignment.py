#%%
import os
#from concurrent.futures import ProcessPoolExecutor
#%%
######################run directly from windows
#def parallel_init(metadata_var, records_var):
#    global metadata
#    metadata = metadata_var
#    global records 
#    records = records_var
def extract_subfamily_coord(subfamily,
                            save_to_file = True, 
                            output_filepath = None, 
                            species_reference=None):
    import pandas as pd
    if species_reference is None:
        species_reference = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/combined_age_div/combined_age_and_div.txt'
    age_div_table = pd.read_csv(species_reference, sep='\t')
    try:
        subfam_table=age_div_table[age_div_table.repName.isin(subfamily)]
    except TypeError:
        print('expecting a list of subfamily')
    subfam_table = subfam_table[subfam_table.genoName.isin(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'])]
    subfam_coord = subfam_table[['genoName','genoStart','genoEnd','strand']].copy()
    subfam_coord.columns = ['chrom','start','end','strand']
    subfam_coord['id'] = subfam_table.repName + '_' + subfam_table.internal_id.astype(str)
    if save_to_file == True:
        if output_filepath is None:
            output_filepath = f'{os.path.dirname(os.path.abspath(__file__))}/{subfamily}.txt'
        subfam_coord.to_csv(output_filepath, sep='\t', index= False)
    else:
        return subfam_coord

def extract_sequence(genoname, genostart, genoend,strand, records):
    #seqname = '::'.join([meta_id,genoname,str(genostart),str(genoend),strand])
    chromosome_extract=records[genoname]
    if strand == '+':
        seq_string = str(chromosome_extract[genostart:genoend].seq)
    else:
        seq_string = str(chromosome_extract[genostart:genoend].seq.reverse_complement())
    return seq_string
    #seq_record = SeqRecord(Seq(''.join(seq_string)),seqname , '', '')
    #return seq_record




def sequence_io(metadata,source_fasta,output_filepath = None, save_to_file=False, custom_id = False):
    import pandas as pd
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    if isinstance(metadata, str):
        if (os.path.isfile(metadata) == True):
            metadata_local = pd.read_csv(metadata, sep='\s+')
        else:
            print('metadata file not found')
    else:
        metadata_local = metadata

    if output_filepath is None:
        if isinstance(metadata, str):
            output_filepath =f"{'/'.join(str.split(metadata, sep ='/')[:-1])}.fasta"
        else:
            output_filepath = f'{os.path.dirname(os.path.abspath(__file__))}/sequence.fasta'
    else:
        output_filepath = output_filepath
    
    records = SeqIO.to_dict(SeqIO.parse(open(source_fasta), 'fasta'))
    seq_records = []
    
    if custom_id == False:
        meta_id = [f'sample_n{index}' for index in metadata_local.index.astype(str)]
        metadata_local['meta_id'] = meta_id
    else:
        metadata_local['meta_id'] = metadata_local.iloc[:,4]
        meta_id = metadata_local.meta_id.unique()
    for uniq_meta_id in meta_id:
        metadata_by_id = metadata_local[metadata_local.meta_id == uniq_meta_id]
        seq_strings = []
        for idx, row in metadata_by_id.iterrows():
            #print(row)       
            chrom = row.chrom
            start = row.start
            end = row.end
            strand = row.strand
            print(chrom, start, end, strand, uniq_meta_id)
            seq_string = extract_sequence(chrom, start, end, strand, records)
            seq_strings.append(seq_string)
        if strand == '-':
            seq_strings.reverse()
        seqname = '::'.join([uniq_meta_id,chrom,str(min(metadata_by_id.iloc[:,1])),str(max(metadata_by_id.iloc[:,2])),strand])
        seq_record = SeqRecord(Seq(''.join(seq_strings)),seqname , '', '')
        seq_records.append(seq_record)

    if save_to_file == True:
        with open(output_filepath, "w") as output_handle:
            SeqIO.write(seq_records, output_handle, "fasta")
    else:
        return seq_records

def parsed_array_to_sequence(parsed_array,prefix, output_filepath):
    import numpy as np
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO

    # Define the mapping
    base_mapping = {0: 'N', 1: 'A', 2: 'C', 3: 'T',4:'G'}

    # Convert each row to DNA sequence
    dna_sequences = []
    for idx, row in enumerate(parsed_array):
        seq_name = f'{prefix}_{idx}'
        dna_sequence = ''.join(base_mapping[int(val)] for val in row)
        dna_sequences.append(SeqRecord(Seq(dna_sequence),seq_name , '', ''))

    # Write the sequences to a FASTA file
    with open(output_filepath, "w") as output_handle:
        SeqIO.write(dna_sequences, output_handle, "fasta")

def mafft_align(input_filepath, nthread = None, nthreadtb = None, nthreadit =None, output_filepath = None, mafft_arg = 
                  #--quiet --memsave 
                  '--treeout --reorder '):
    import subprocess
    mafft_command = 'mafft '
    if nthread is not None:
        nthread_arg = '--thread '+ str(nthread) + ' ' 
        mafft_command += nthread_arg
    if nthreadtb is not None:
        nthreadb_arg = '--threadb '+ str(nthread) + ' ' 
        mafft_command += nthreadb_arg
    if nthreadit is not None:
        nthreadit_arg = '--threadit '+ str(nthread) + ' ' 
        mafft_command += nthreadit_arg
    if mafft_arg is not None:
        mafft_command += mafft_arg
    mafft_command += input_filepath
    if output_filepath is None:
        output_filepath = f'{input_filepath}.aligned' 
    try:
        output = subprocess.check_output(mafft_command, shell = True)
        with open(output_filepath, 'wb') as file:
            file.write(output)
    except subprocess.CalledProcessError as e:
        print(e.cmd)
        print(e.output)
        print(e.returncode)
        print("MAFFT error")
        raise
#%%
