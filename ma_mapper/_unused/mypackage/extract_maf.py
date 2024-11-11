#%%
import numpy as np
import pandas as pd
import cyvcf2
import os
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat 
from Bio.AlignIO import MafIO
if sys.version_info >= (3, 8, 0):
    from typing import Literal
else:
    from typing_extensions import Literal
#%%
_AGEARG = Literal[None,'extract','calibrate']
_COUNTARG = Literal['human_ref','coverage','common', 'common_raw', 'total_raw']
def extract_maf(maf_file, chrom, start, end, strand,target_species = "Homo_sapiens", count_arg:_COUNTARG = 'human_ref',age=None, age_arg:_AGEARG = None, age_table_file = None):
    #print(target_species,chrom, start, end, strand)
    if age_arg is not None:
        if isinstance(age_table_file, str):
            if (os.path.isfile(age_table_file) == True):
                species = pd.read_table(age_table_file, sep='\t')
            else:
                print('metadata file not found')
        else:
            species = age_table_file
    maf_id = target_species+'.'+chrom

    index_maf = MafIO.MafIndex(maf_file+".mafindex", maf_file, maf_id) 
    if strand =='-':
        n_strand = -1
        results =index_maf.get_spliced(start,end,n_strand)
    else:
        n_strand = 1
        results =index_maf.get_spliced(start,end,n_strand)
    collector = {}
    for seqrec in results:
        if seqrec.id.split('.')[0] == target_species:  
            collector[seqrec.id.split('.')[0]]= np.array(list(seqrec.seq.lower()))
            age_measure = 0
            for seqrec in results:
                if seqrec.id.split('.')[0] != target_species:
                    if age_arg is None:
                        collector[seqrec.id.split('.')[0]]= np.array(list(seqrec.seq.lower()))
                    elif age_arg == 'calibrate':
                        species_age = species[species.meta_name == seqrec.id.split('.')[0]]['Estimated Time (MYA)'].to_list()[0]
                        if age >= species_age:
                            collector[seqrec.id.split('.')[0]]= np.array(list(seqrec.seq.lower()))
                    elif age_arg == 'extract':
                        species_age = species[species.meta_name == seqrec.id.split('.')[0]]['Estimated Time (MYA)'].to_list()[0]
                        if species_age >= age_measure:
                            age_measure = species_age
    if age_arg == 'extract':
        return age_measure
    else:
        array_transposed=np.array(list(collector.values())).transpose()
        output_array=[]
        for ref_pos in array_transposed:
            (unique, counts) = np.unique(ref_pos, return_counts=True)
            frequencies = dict(zip(unique, counts))
            ref_allele=ref_pos[0]
            #frequencies.pop('-', None) #->count deletion/insertion
            total = sum(frequencies.values())
            if count_arg == 'human_ref':
                alt_count = total - frequencies[ref_allele]
                alt_freq=alt_count/total
                output_array.append(alt_freq)
            elif count_arg == 'coverage':
                output_array.append(total)
            elif count_arg == 'common':
                common_allele = max(frequencies, key=frequencies.get)
                common_count = frequencies[common_allele]
                common_freq = common_count/total
                output_array.append(common_freq)
            elif count_arg == 'common_raw':
                common_allele = max(frequencies, key=frequencies.get)
                common_count = frequencies[common_allele]
                output_array.append(common_count)
            elif count_arg == 'total_raw':
                output_array.append(total)
        return output_array
#%%
def fetch_maf(metadata_input, maf_input,output_dir = None, separated_maf = False, target_species = 'Homo_sapiens', count_arg = 'human_ref', save_to_file = False, custom_id = False, age_arg:_AGEARG = None, age_table_file = None, age_list=None):
    if isinstance(metadata_input, str):
        if os.path.isfile(metadata_input):
            metadata = pd.read_csv(metadata_input, sep='\t')
        else:
            print('metadata file not found')
    else:
        metadata = metadata_input
    if output_dir is None:
        if isinstance(metadata_input, str):
            output_dir = '/'.join(str.split(metadata_input, sep ='/')[:-1])
        else: 
            output_dir = os.path.dirname(os.path.abspath(__file__))

    print('extract from maf target: '+ maf_input)
    if custom_id == False:
        meta_id = 'sample_n'+ metadata.index.astype(str)
        metadata['meta_id'] = meta_id
    else:
        metadata['meta_id'] = metadata.iloc[:,4]
        meta_id = metadata.iloc[:,4].unique()
    maf_call_list = []
    chrom_list = []
    start_list = []
    end_list = []
    strand_list = []
    for uniq_meta_id in meta_id:
        metadata_by_id = metadata[metadata.meta_id == uniq_meta_id]
        chrom = metadata_by_id.iloc[:,0].unique()[0]
        chrom_list.append(chrom)
        start_list.append(metadata_by_id.iloc[:,1].to_list())
        end_list.append(metadata_by_id.iloc[:,2].to_list())
        strand_list.append(metadata_by_id.iloc[:,3].unique()[0])
        if separated_maf == True:
            maf_file = maf_input + '.' + chrom
        else:
            maf_file = maf_input
        maf_call_list.append(maf_file)
    if age_arg is None:
        print('normal maf extraction without calibration')
        with ProcessPoolExecutor(max_workers=40) as executor:
            results  = executor.map(extract_maf, maf_call_list, chrom_list, start_list, end_list, strand_list, repeat(target_species) ,repeat(count_arg))
    elif age_arg == 'extract':
        print('extract age depth for further use')
        with ProcessPoolExecutor(max_workers=40) as executor:
            results  = executor.map(extract_maf, maf_call_list, chrom_list, start_list, end_list, strand_list, repeat(target_species) ,repeat(count_arg), repeat(age_list), repeat(age_arg), repeat(age_table_file))
    elif age_arg == 'calibrate':
        if age_list is None:
            print('calibrate age: no age_list provided, fetching from metadata')
            age_list = []
            for uniq_meta_id in meta_id:
                age_list.append(metadata_by_id.te_age.unique()[0])
        else:
            print('calibrate age: using preexisting age_list')
            age_list = age_list
        with ProcessPoolExecutor(max_workers=40) as executor:
            results  = executor.map(extract_maf, maf_call_list, chrom_list, start_list, end_list, strand_list, repeat(target_species) ,repeat(count_arg), age_list, repeat(age_arg), repeat(age_table_file))
    maf_out = []
    for result in results:
        maf_out.append(result)
    if save_to_file == True:
        import compress_pickle
        output_filepath = output_dir+'/maf_out.lzma'
        compress_pickle.dump(maf_out, output_filepath, compression="lzma")
        print('done, saving maf_out at: '+output_filepath)
    else:
        print('done, returning maf_out as object')
        return maf_out
#%%