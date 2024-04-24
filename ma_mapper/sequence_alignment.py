#%%
import os
global species_reference
species_reference = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/combined_age_div/combined_age_and_div.txt'
global main_chr
main_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
#from concurrent.futures import ProcessPoolExecutor
#%%
######################run directly from windows
#def parallel_init(metadata_var, records_var):
#    global metadata
#    metadata = metadata_var
#    global records 
#    records = records_var
def extract_subfamily_coord(subfamily,save_to_file = True, output_filepath = None):
    import pandas as pd
    age_div_table = pd.read_csv(species_reference, sep='\t')
    try:
        subfam_table=age_div_table[age_div_table.repName.isin(subfamily)]
    except TypeError:
        print('expecting a list of subfamily')
    subfam_table = subfam_table[subfam_table.genoName.isin(main_chr)]
    subfam_coord = subfam_table[['genoName','genoStart','genoEnd','strand']].copy()
    subfam_coord['id'] = subfam_table.repName + '_' + subfam_table.internal_id.astype(str)
    if save_to_file == True:
        if output_filepath is None:
            output_filepath = os.path.dirname(os.path.abspath(__file__))
        subfam_coord.to_csv(output_filepath, sep='\t', index= False)
    else:
        return subfam_coord

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

def extract_sequence(genoname, genostart, genoend, strand):
    chromosome_extract=records[genoname]
    if strand == '+':
        seq_string = str(chromosome_extract[genostart:genoend].seq)
    else:
        seq_string = str(chromosome_extract[genostart:genoend].seq.reverse_complement())
    return seq_string

def fetch_sequence(metadata_input,source_fasta,output_filepath = None, save_to_file=False, custom_id = False):
    import pandas as pd
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    if isinstance(metadata_input, str):
        if (os.path.isfile(metadata_input) == True):
            metadata = pd.read_csv(metadata_input, sep='\s+')
        else:
            print('metadata file not found')
    else:
        metadata = metadata_input

    if output_filepath is None:
        if isinstance(metadata_input, str):
            output_filepath = '/'.join(str.split(metadata_input, sep ='/')[:-1]) 
        else:
            output_filepath = os.path.dirname(os.path.abspath(__file__))
    else:
        output_filepath = '/'.join(str.split(output_filepath, sep ='/')[:-1]) 
    global records
    records = SeqIO.to_dict(SeqIO.parse(open(source_fasta), 'fasta'))
    seq_records = []
    
    if custom_id == False:
        meta_id = 'sample_n'+ metadata.index.astype(str)
        metadata['meta_id'] = meta_id
    else:
        metadata['meta_id'] = metadata.iloc[:,4]
        meta_id = metadata.iloc[:,4].unique()

    for uniq_meta_id in meta_id:
        metadata_by_id = metadata[metadata.meta_id == uniq_meta_id]
        seq_strings = []
        for idx, row in metadata_by_id.iterrows():       
            chrom = row.chrom
            start = row.start
            end = row.end
            strand = row.strand
            seq_string = extract_sequence(chrom, start, end, strand)
            seq_strings.append(seq_string)
        seqname = '::'.join([uniq_meta_id,chrom,str(min(metadata_by_id.iloc[:,1])),str(max(metadata_by_id.iloc[:,2])),strand])
        seq_record = SeqRecord(Seq(''.join(seq_strings)),seqname , '', '')
        seq_records.append(seq_record)

    if save_to_file == True:
        with open(output_filepath, "w") as output_handle:
            SeqIO.write(seq_records, output_handle, "fasta")
    else:
        return seq_records
    
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
        output_filepath = input_filepath+'.aligned' 
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
