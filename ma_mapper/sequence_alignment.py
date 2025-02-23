#%%
import os
import pandas as pd
#from concurrent.futures import ProcessPoolExecutor
#%%
######################run directly from windows
#def parallel_init(coordinate_var, records_var):
#    global coordinate_table
#    coordinate_table = coordinate_var
#    global records 
#    records = records_var
def extract_subfamily_coord_from_combined_table(subfamily,
                            save_to_file = True, 
                            output_filepath = None, 
                            species_reference=None):
    
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
#%%

def extract_coord_from_repeatmasker_table(
        subfamily: str|list,
        repeatmasker_table: str|pd.DataFrame,
        internal_id_table: str|list|None=None,
        save_to_file:bool = False,
        output_filepath: str|None = None,
        ):
    
    if isinstance(subfamily, str):
        subfamilies = [subfamily]
    else:
        subfamilies = subfamily
    coord_df = []
    for idx,subfamily in enumerate(subfamilies):
        if isinstance(repeatmasker_table, str):
            repeatmasker_table = pd.read_csv(repeatmasker_table, sep='\t', index_col=0, header=0)

        if internal_id_table is not None:
            if isinstance(internal_id_table, list):
                internal_id_tbl = internal_id_table[idx]
            else:
                internal_id_tbl = internal_id_table
            if isinstance(internal_id_tbl, str):
                if os.path.isfile(internal_id_tbl):
                    internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
                else:
                    raise FileNotFoundError(f"Internal ID table {internal_id_tbl} not found")

            elif isinstance(internal_id_tbl, pd.DataFrame):
                internal_id_df = internal_id_tbl
            subfam_w_internal_id=pd.merge(repeatmasker_table, internal_id_df, left_index = True, right_on = 'rmsk_index')
            subfam_coords = subfam_w_internal_id[['genoName','genoStart','genoEnd','internal_id']].copy()
            subfam_coords['score'] = 10
            subfam_coords['strand'] = subfam_w_internal_id.strand
        else:
            subfam_table=repeatmasker_table[repeatmasker_table.repName == subfamily].copy().reset_index()
            subfam_coords = subfam_table[['genoName','genoStart','genoEnd']].copy()
            subfam_coords['internal_id'] =subfam_table.index.map(lambda x: f'{subfamily}_{x}')
            subfam_coords['score'] = 10
            subfam_coords['strand'] = subfam_table.strand
        coord_df.append(subfam_coords)
    #     == True:
    coords = pd.concat(coord_df).reset_index(drop=True)
    if save_to_file== True:
        if output_filepath is None:
            output_filepath = f'{os.path.dirname(os.path.abspath(__file__))}/{subfamily}.txt'
        coords.to_csv(output_filepath, sep='\t', index= False, header=False)
    return coords
#%%
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




def sequence_io(coordinate_table,
                source_fasta,
                output_filepath = None, 
                save_to_file=False, 
                custom_id = None, 
                ):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    if isinstance(coordinate_table, str):
        if (os.path.isfile(coordinate_table) == True):
            coordinate_local = pd.read_csv(coordinate_table, sep='\t', header=None)
        else:
            print('coordinate_table file not found')
    else:
        coordinate_local = coordinate_table

    if output_filepath is None:
        if isinstance(coordinate_table, str):
            output_filepath =f"{'/'.join(str.split(coordinate_table, sep ='/')[:-1])}.fasta"
        else:
            output_filepath = f'{os.path.dirname(os.path.abspath(__file__))}/sequence.fasta'
    else:
        output_filepath = output_filepath
    
    if custom_id:
        meta_id = [f'{custom_id}_{index}' for index in coordinate_local.index.astype(str)]
        coordinate_local['meta_id'] = meta_id
    else:
        coordinate_local['meta_id'] = coordinate_local.iloc[:,3]
        meta_id = coordinate_local.meta_id.unique()

    records = SeqIO.to_dict(SeqIO.parse(open(source_fasta), 'fasta'))
    seq_records = []
    for uniq_meta_id in meta_id:
        coordinate_by_id = coordinate_local[coordinate_local.meta_id == uniq_meta_id]
        seq_strings = []
        for idx, row in coordinate_by_id.iterrows():
            #print(row)       
            chrom = row.iloc[0]
            start = row.iloc[1]
            end = row.iloc[2]
            strand = row.iloc[5]
            #print(chrom, start, end, strand, uniq_meta_id)
            seq_string = extract_sequence(chrom, start, end, strand, records)
            seq_strings.append(seq_string)
        if strand == '-':
            seq_strings.reverse()
        seqname = f'{uniq_meta_id}::{chrom}:{min(coordinate_by_id.iloc[:,1])}-{max(coordinate_by_id.iloc[:,2])}({strand})'
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
