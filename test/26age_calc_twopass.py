#%%
import pandas as pd
import logging
from datetime import datetime
import os
import time
import itertools
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import sys
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/stable')
import e_val_calc
import config_hg38_repeatlib as config
from Bio.AlignIO import MafIO
#%%
def age_filter(internal_id,e_cutoff = 1e-3):
    temp_tbl = internal_id_tbl[internal_id_tbl.internal_id == internal_id]
    cutoff_tbl = temp_tbl[temp_tbl.E_value.astype('float64') <= e_cutoff]
    temp_age=cutoff_tbl.age.max()
    return temp_age
#%%
def batch_age_calc(subfamily):
    subfamily_filename = subfamily.replace('/','_')
    output_filepath = output_folder+'/'+subfamily_filename+'.txt'
    operation_log_path = output_folder+'/op_log.log'
    logging.root.handlers = []
    logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    handlers = [
                        logging.FileHandler(operation_log_path, mode = 'a'),
                        logging.StreamHandler()
                        ]
                    )
    logging.info('process: '+subfamily)
    logging.debug('process: '+output_filepath)
    age_list = list()
    if (os.path.isfile(output_filepath) == False):
        start_time = time.time()
        internal_id_list = internal_id_tbl.internal_id.unique()
        with ProcessPoolExecutor(max_workers=40) as executor:
            results = list(tqdm(executor.map(age_filter, internal_id_list ,itertools.repeat(e_cutoff)), total=len(internal_id_list)))
            #results = executor.map(age_filter, internal_id_list,itertools.repeat(internal_id_tbl) ,itertools.repeat(e_cutoff))
        for result in results:
            age_list.append(result)
        dict_prep = {'internal_id': internal_id_list, 'te_age':age_list,}
        output_table=pd.DataFrame(dict_prep)
        output_table.to_csv(output_filepath, sep='\t', index=False)
        logging.info('done, saving the table at: '+output_filepath)
        logging.debug("--- %s seconds ---\n" % (time.time() - start_time))
    else:
        logging.info(subfamily+' already done')
#%%
def main():
    global e_cutoff, input_folder, output_folder
    e_cutoff = 1e-3
    input_folder = config.e_value_folder
    output_folder = config.te_age_folder
    if os.path.isdir(output_folder) == False:
        os.mkdir(output_folder) 
    for subfamily in config.subfamily_list:
        subfamily_filename = subfamily.replace('/','_')
        input_filepath = input_folder+'/'+subfamily_filename+'.txt'
        global internal_id_tbl
        internal_id_tbl = pd.read_csv(input_filepath, sep='\t', low_memory=False)
        batch_age_calc(subfamily)
#%%
e_cutoff = 1e-3
subfamily = 'MER11C'
input_folder = config.internal_id_folder
subfamily_filename = subfamily.replace('/','_')
input_filepath = input_folder+'/'+subfamily_filename+'.txt'
internal_id_table = pd.read_csv(input_filepath, sep='\t', low_memory=False)
internal_id_list = internal_id_table.internal_id.unique()

internal_id = 320
#%%
rmsk_index=internal_id_table[internal_id_table.internal_id == internal_id].rmsk_index
filtered_table = config.filtered_table
table_by_id= filtered_table[filtered_table.index.isin(rmsk_index)]

chrom=table_by_id.genoName.unique()[0]
strand=table_by_id.strand.unique()[0]
genostart_list=table_by_id.genoStart.to_list()
genoend_list=table_by_id.genoEnd.to_list()
full_length = [min(genostart_list)-5000,max(genoend_list)+5000]
if strand=='-':
    strand = -1
else:
    strand = 1
target_species = config.target_species
species_table = config.species_table
if target_species == 'Homo_sapiens':
    mafpath = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/241genomes/241-mammalian-2020v2b.maf.'+ chrom
elif target_species == 'Mus_musculus':
    mafpath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/mouse241/Mus_musculus.maf.'+ chrom
target_chrom = target_species+'.'+chrom
index_maf = MafIO.MafIndex(mafpath+".mafindex", mafpath, target_chrom)
#%%
spliced_maf_full =index_maf.get_spliced([full_length[0]],[full_length[1]],strand)
collector = []
for seqrec in spliced_maf_full:
    if seqrec.id.split('.',maxsplit=1)[0] == target_species:
        target_seq = seqrec.seq
        for seqrec in spliced_maf_full:
            meta_name = seqrec.id.split('.',maxsplit=1)[0]
            species_id = species_table[species_table['meta_name']==meta_name].index[0]
            chrom_name = seqrec.id.split('.',maxsplit=1)[1]
            #first pass - check TE E-value
            count_result = e_val_calc.affine_count_simple(target_seq[5000:-5000], seqrec.seq[5000:-5000], matchWeight = 1, mismatchWeight = 1,            gapWeight = 4, extendedgapWeight = 1,)
            blast_score=count_result[0]; matched = count_result[1]; gapped =count_result[2]; gap_count = count_result[3]
            n = species_table[species_table['meta_name']==seqrec.id.split('.',maxsplit=1)[0]]['ungapped_length'].values[0]
            m = len(seqrec.seq)
            if target_species == 'Homo_sapiens':
                age=species_table[species_table['meta_name']==seqrec.id.split('.',maxsplit=1)[0]]['Estimated Time (MYA)'].values[0]
            elif target_species == 'Mus_musculus':
                age=species_table[species_table['meta_name']==seqrec.id.split('.',maxsplit=1)[0]]['Estimated Time mouse_ref (MYA)'].values[0]
            if matched != 0:
                p_value, E = e_val_calc.BLAST_StoP(blast_score, m,n ,lambda_=1.08,K=0.28, H=0.54, alpha=2.0, beta=-2,gapped = gapped)
                if E <= e_cutoff:
                    #second pass - check flanking regions E-value
                    count_result_front = e_val_calc.affine_count_simple(target_seq[:5000], seqrec.seq[:5000], matchWeight = 1, mismatchWeight = 1,gapWeight = 4, extendedgapWeight = 1,)
                    blast_score_front=count_result_front[0]; matched_front = count_result_front[1]; gapped_front =count_result_front[2]; gap_count_front = count_result_front[3]
                    p_value_front, E_front = e_val_calc.BLAST_StoP(blast_score_front, m,n ,lambda_=1.08,K=0.28, H=0.54, alpha=2.0, beta=-2,gapped = gap_count_front)
                    count_result_back = e_val_calc.affine_count_simple(target_seq[-5000:], seqrec.seq[-5000:], matchWeight = 1, mismatchWeight = 1,gapWeight = 4, extendedgapWeight = 1,)
                    blast_score_back=count_result_back[0]; matched_back = count_result_back[1]; gapped_back =count_result_back[2]; gap_count_back = count_result_back[3]
                    p_value_back, E_back = e_val_calc.BLAST_StoP(blast_score_back, m,n ,lambda_=1.08,K=0.28, H=0.54, alpha=2.0, beta=-2,gapped = gap_count_back)
                    if (E_front <= e_cutoff or E_back <= e_cutoff):
                        iden=round((matched/len(target_seq[5000:-5000])*100), 2)
                        pgap=round((gap_count/len(target_seq[5000:-5000])*100), 2)
                        E_score ='{:0.2e}'.format(E)
                        iden_front=round((matched_front/5000*100), 2)
                        pgap_front=round((gap_count_front/5000*100), 2)
                        E_score_front ='{:0.2e}'.format(E_front)
                        iden_back=round((matched_back/5000*100), 2)
                        pgap_back=round((gap_count_back/5000*100), 2)
                        E_score_back ='{:0.2e}'.format(E_back)
                        collector.append([internal_id,meta_name,chrom_name,age, iden, pgap, blast_score, E_score,[iden_front,iden_back],[pgap_front,pgap_back],[E_score_front,E_score_back]])
# post-hoc flank region check for 'human specific'
if len(collector) == 1:
    post_hoc_collector = []
    match_count = 0; total_count = 0; front_only_count = 0; back_only_count = 0; both_count=0; nonmatch_count = 0;
    for seqrec in spliced_maf_full:
        if seqrec.id.split('.',maxsplit=1)[0] == target_species:
            target_seq = seqrec.seq
            for seqrec in spliced_maf_full:
                #second pass - check flanking regions E-value
                count_result_front = e_val_calc.affine_count_simple(target_seq[:5000], seqrec.seq[:5000], matchWeight = 1, mismatchWeight = 1,gapWeight = 4, extendedgapWeight = 1,)
                blast_score_front=count_result_front[0]; matched_front = count_result_front[1]; gapped_front =count_result_front[2]; gap_count_front = count_result_front[3]
                p_value_front, E_front = e_val_calc.BLAST_StoP(blast_score_front, m,n ,lambda_=1.08,K=0.28, H=0.54, alpha=2.0, beta=-2,gapped = gap_count_front)
                count_result_back = e_val_calc.affine_count_simple(target_seq[-5000:], seqrec.seq[-5000:], matchWeight = 1, mismatchWeight = 1,gapWeight = 4, extendedgapWeight = 1,)
                blast_score_back=count_result_back[0]; matched_back = count_result_back[1]; gapped_back =count_result_back[2]; gap_count_back = count_result_back[3]
                p_value_back, E_back = e_val_calc.BLAST_StoP(blast_score_back, m,n ,lambda_=1.08,K=0.28, H=0.54, alpha=2.0, beta=-2,gapped = gap_count_back)
                total_count +=1
                if (E_front <= e_cutoff or E_back <= e_cutoff):
                    match_count += 1
                    if E_front <= e_cutoff and E_back >= e_cutoff:
                        front_only_count += 1
                    elif E_front >= e_cutoff and E_back <= e_cutoff:
                        back_only_count +=1
                    elif E_front <= e_cutoff and E_back <= e_cutoff:
                        both_count += 1
                else:
                    nonmatch_count += 1
    post_hoc_collector = collector[0][:8]
    post_hoc_collector.extend([[match_count,total_count],[front_only_count,back_only_count],[both_count,nonmatch_count]])
#%% convert to df
pd.DataFrame(collector, columns= ['internal_id','species','chr_code','divergence','%iden','%gap','BLAST_score','E_value','%iden_flanks','%gap_flank','E_val_flanks']).sort_values('divergence',ascending =True)
# %%
pd.DataFrame([post_hoc_collector], columns= ['internal_id','species','chr_code','divergence','%iden','%gap','BLAST_score','E_value','%iden_flanks','%gap_flank','E_val_flanks']).sort_values('divergence',ascending =True)
# %%
