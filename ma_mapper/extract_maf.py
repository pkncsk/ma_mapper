#%%
import numpy as np
import pandas as pd
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat 
from Bio.AlignIO import MafIO
from . import logger
if sys.version_info >= (3, 8, 0):
    from typing import Literal, Tuple, List
else:
    from typing_extensions import Literal, Tuple, List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter, defaultdict
import time

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

MafIO.MafIndex.get_spliced = get_spliced_mod
#%%
_AGEARG = Literal[None,'calibrate']
_COUNTARG = Literal['ref_freq','coverage','common_freq', 'common_raw','common_freq_nogap','a','t','c','g','base_count','base_freq','base_freq_nogap','raw']
def extract_maf(name:str,
                maf_file:str, 
                chrom:str, 
                start:int, 
                end:int, 
                strand:str,
                target_species:str = 'Homo_sapiens', 
                count_arg:_COUNTARG = 'common',
                species_list:list=None,
                e_value_df:pd.DataFrame = None, 
                internal_id_df:pd.DataFrame = None,
                ):
    iupac_codes = {
        'A': {'A'},
        'C': {'C'},
        'G': {'G'},
        'T': {'T'},
        'R': {'A', 'G'},
        'Y': {'C', 'T'},
        'S': {'G', 'C'},
        'W': {'A', 'T'},
        'K': {'G', 'T'},
        'M': {'A', 'C'},
        'B': {'C', 'G', 'T'},
        'D': {'A', 'G', 'T'},
        'H': {'A', 'C', 'T'},
        'V': {'A', 'C', 'G'},
        'N': {'A', 'C', 'G', 'T'},
        '-': {'-'}
    }

    def count_bases_with_ambiguity(sequences):
        base_counts = Counter()
        for sequence in sequences:
            for base in sequence:
                possible_bases = iupac_codes.get(base.upper(), {base.upper()})
                # Distribute counts among all possible bases
                for possible_base in possible_bases:
                    base_counts[possible_base] += 1 / len(possible_bases)
        
        return base_counts

    maf_id = f'{target_species}.{chrom}'
    
    index_maf = MafIO.MafIndex(f'{maf_file}.mafindex', maf_file, maf_id) 
    n_strand = -1 if strand == '-' else 1
    results =index_maf.get_spliced(start,end,n_strand)
    if e_value_df is None and internal_id_df is None and species_list is None:
        collector = results
    elif species_list is not None:
        pattern = '|'.join(species_list)  
        collector = results[results['seqid'].str.contains(pattern)]
    else:
        if internal_id_df is not None:
            _internal_id = internal_id_df[internal_id_df['meta_id'] == internal_id]['internal_id'].values[0]
        else:
            _internal_id = name
        e_value_df['seqid'] = e_value_df.species + '.' + e_value_df.chr_code
        e_value_internal_id=e_value_df[e_value_df['internal_id'] == _internal_id]
        
        if e_value_internal_id.shape[0] <= 1:
            collector = results[results.seqid.str.contains(target_species)]
        else:
            collector = results[results.seqid.isin(e_value_internal_id.seqid.values)]
            
    if count_arg == 'raw':
        sequence_length = len(results.iloc[0]['seq'])
        list_of_dfs = []
        for i in range(sequence_length):
            # Extract the base at the i-th position for each row
            temp_df = results[['seqid']].copy()  # Start with the seqid column
            temp_df['seq'] = results['seq'].apply(lambda x: x[i])  # Add the base at the i-th position
            list_of_dfs.append(temp_df)
        return list_of_dfs
    #return collector
    #logger.info('pass1')
    try:

        ref_alleles = np.char.upper(results[results.seqid.str.contains(target_species)]['seq'].to_list())[0]
    except IndexError:
        print(f'IndexError:{name}\t{maf_file}\t{maf_id}\t{chrom}\t{start}\t{end}\t{strand}')
       
    
    
    array_transposed=np.array(collector['seq'].to_list()).transpose()
    output_array=[]
    for idx, pos_array in enumerate(array_transposed):
        frequencies =count_bases_with_ambiguity(np.char.upper(pos_array))
        ref_allele=ref_alleles[idx]
        total = sum(frequencies.values())
        if count_arg == 'ref_freq':
            alt_count = total - frequencies[ref_allele]
            alt_freq=alt_count/total
            if total == 1:
                alt_freq = np.nan
            output_array.append(alt_freq)
        elif count_arg == 'coverage':
            output_array.append(total)
        elif count_arg in ['common_freq', 'common_raw']:
            common_allele = max(frequencies, key=frequencies.get)
            common_count = frequencies.get(common_allele)
            if count_arg == 'common_freq':
                common_freq = common_count/total
                output_array.append(common_freq)
            else:
                output_array.append(common_count)
        elif count_arg == 'common_freq_nogap':
            frequencies.pop('-', None)
            try:
                common_allele = max(frequencies, key=frequencies.get)
            except ValueError:
                print(f'ValueError:{name}\t{maf_file}\t{maf_id}\t{chrom}\t{start}\t{end}\t{strand}')
            common_count = frequencies.get(common_allele)
            common_freq = common_count/total
            output_array.append(common_freq)
        elif count_arg == 'base_count':
            output_array.append(frequencies)
        elif count_arg in ['base_freq','base_freq_nogap']:
            if count_arg == 'base_freq_nogap':
                frequencies.pop('-', None)
            
            count_freq = {base: count / total for base, count in frequencies.items()}
            output_array.append(count_freq)
    if isinstance(output_array, str):
        print('debug',name,output_array)
    return output_array
#%%
def maf_io(coordinate_table: pd.DataFrame|str, 
           maf:str,
           separated_maf:bool = False, 
           target_species:str = 'Homo_sapiens', 
           count_arg:_COUNTARG = 'common', 
           save_to_file:bool = False, 
           generate_new_id:bool = False, 
           species_list:list=None,
           e_value_table: str|pd.DataFrame|None=None, 
           internal_id_table: str|pd.DataFrame|None=None, **kwargs):
    if isinstance(coordinate_table, str):
        if os.path.isfile(coordinate_table):
            coordinate_local = pd.read_csv(coordinate_table, sep='\t', header=None)
        else:
            logger.error('coordinate_table file not found')
    else:
        coordinate_local = coordinate_table
    

    logger.info(f'extract from maf target: {maf}')
    if generate_new_id == True:
        meta_id = [f'entry_{index}' for index in coordinate_local.index.astype(str)]
        coordinate_local['meta_id'] = meta_id
    else:
        coordinate_local['meta_id'] = coordinate_local.iloc[:,3]
        meta_id = coordinate_local.meta_id.unique()
    print(meta_id[0])
    grouped = coordinate_local.groupby('meta_id', sort=False)
    chrom_list = grouped.apply(lambda x: x.iloc[:,0].unique()[0]).tolist()
    start_list = grouped.apply(lambda x: x.iloc[:,1].tolist()).tolist()
    end_list = grouped.apply(lambda x: x.iloc[:,2].tolist()).tolist()
    strand_list = grouped.apply(lambda x: x.iloc[:,5].unique()[0]).tolist()
    maf_call_list = []
    for chrom in chrom_list:
        if separated_maf == True:
            maf_file = f'{maf}.{chrom}'
        else:
            maf_file = maf
        maf_call_list.append(maf_file)
    e_value_df = None
    if isinstance(e_value_table, str):
        e_value_df = pd.read_csv(e_value_table, sep='\t')
    else:
        e_value_df = e_value_table
    internal_id_df = None
    if isinstance(internal_id_table, str):
        internal_id_df = pd.read_csv(internal_id_table, sep='\t')
    else:
        internal_id_df = internal_id_table

    with ProcessPoolExecutor(max_workers=40) as executor:
        results = executor.map(extract_maf, meta_id, maf_call_list, chrom_list, start_list, end_list, strand_list, repeat(target_species), repeat(count_arg), repeat(e_value_df), repeat(internal_id_df), repeat(species_list))

    maf_out = []
    for result in results:
        maf_out.append(result)
    if save_to_file == True:
        if isinstance(save_to_file, str):
            output_filepath = save_to_file
        else:
            if isinstance(coordinate_table, str):
                output_dir = '/'.join(str.split(coordinate_table, sep ='/')[:-1])
            else: 
                output_dir = os.path.dirname(os.path.abspath(__file__))
            output_filepath = f'{output_dir}/maf_output.p'
        import compress_pickle
        compress_pickle.dump(maf_out, output_filepath, compression="lzma")
        logger.info('done, saving maf_out at: ', output_filepath)
    else:
        logger.info('done, returning maf_out as object')
        return maf_out
#%%
