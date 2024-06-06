# NOTE: 30 August add main(), redirect common used element to config.py
# NOTE: 7 August v1.1 change common length calculation from genoEnd-genoStart to repEnd when repLeft = 0
# NOTE: 4 July completed version of TE position code
#%%
import time
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import sys
import os
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/stable')
import config_hg38_repeatlib as config
#import config_mm39_dfam as config
#import config
import logging
#%%
def assign_te_internal_id(subfamily):
    start_time = time.time()
    filtered_table =config.filtered_table
    operation_log_path = output_folder+'/op_log.log'
    #setup logger
    logging.root.handlers = []
    logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    handlers = [
                        logging.FileHandler(operation_log_path, mode = 'a'),
                        logging.StreamHandler()
                        ]
                    )
    logging.info('assign internal id to te: greedy method')
    subfamily_filename = subfamily.replace('/','_')
    output_filepath = output_folder+'/'+subfamily_filename+'.tsv'
    logging.info('process:'+subfamily)
    if (os.path.isfile(output_filepath) == True):
        logging.info('already done')
    else:
        subfam_table=filtered_table[filtered_table['repName']==subfamily]
        # find common length of TE
        subfam_table=subfam_table.assign(length = subfam_table.genoEnd-subfam_table.genoStart)
        subfam_table=subfam_table.assign(repLength = subfam_table.repEnd-subfam_table.repStart+1)
        #check common length availability
        if subfam_table[subfam_table.repLeft==0].shape[0] != 0:
            common_length = subfam_table[subfam_table.repLeft==0].repEnd.mode()[0]
        else: 
            common_length = (subfam_table.repLeft + subfam_table.repEnd).mode()[0]
        #common_length=subfam_table.length.mode()[0]
        # find unique ID according to repeatmasker 
        # look at subfam table, find entires with 1 fragment (separate them for ID assignment later)
        #subfam_table
        row_count = subfam_table.id.value_counts()
        singlet_id=row_count.index[row_count.lt(2)].to_list()
        singlet_index=subfam_table[subfam_table.id.isin(singlet_id)].index.tolist()
        # look for TE group which has all full length members
        mult_frag_id=row_count.index[row_count.gt(1)].to_list()
        all_complete_frag_id = []
        partial_frag_id = []
        all_complete_frag_index = []
        for idx in mult_frag_id:
            subfam_by_id_tbl = subfam_table[subfam_table.id == idx]
            if all(subfam_row.repLength>common_length * 0.95 for subfam_idx, subfam_row in subfam_by_id_tbl.iterrows()):
                #print(subfam_by_id_tbl.repLength)
                all_complete_frag_id.append(idx)
                all_complete_frag_index.extend(subfam_by_id_tbl.index.tolist())
            else:
                partial_frag_id.append(idx)
        no_mergr_complete_frag_index = []
        no_merge_frag_index = []
        merge_frag_index = []
        log_filepath = output_folder+'/'+subfamily_filename+'.log'
        print(subfamily,'log', sep ='\t', file=open(log_filepath, 'w'))
        for idx in partial_frag_id:
            table_by_id = filtered_table[filtered_table.id == idx]
            table_strand =table_by_id.strand.unique()[0]
            if table_strand =='-':
                table_by_id  = table_by_id.sort_index(ascending=False)
            merge_switch = False
            for row_num, (table_idx,table_row) in enumerate(table_by_id.iterrows()):

                curr_repStart = table_row.repStart
                curr_repEnd = table_row.repEnd
                curr_repLeft = table_row.repLeft
                genoLength = table_row.genoEnd - table_row.genoStart
                repLength = curr_repEnd - curr_repStart
                if table_row.repName == subfamily:
                    if repLength <= common_length * 0.95:
                        # not last row
                        if row_num < table_by_id.shape[0]-1:
                            if not merge_switch:
                                next_repName = table_by_id.iloc[row_num+1].repName
                                if next_repName != subfamily:
                                    merge_switch = False
                                    print(table_row.id, table_row._name, curr_repStart, curr_repEnd,table_row.repLeft,'no_merge', sep ='\t', file=open(log_filepath, 'a'))
                                    no_merge_frag_index.append(table_row._name)
                                else:
                                    next_repStart = table_by_id.iloc[row_num+1].repStart
                                    next_repEnd = table_by_id.iloc[row_num+1].repEnd
                                    overlap_covr = round((curr_repEnd - next_repStart+1)/genoLength*100, 2)
                                    if (overlap_covr < overlap_cov_threshold):
                                        #& (next_repLength <= common_length * 0.95)
                                        merge_switch = True
                                        print(table_row.id, table_row._name, curr_repStart, curr_repEnd,table_row.repLeft,'start', overlap_covr, sep ='\t', file=open(log_filepath, 'a'))
                                        merge_set = [table_row._name]
                                    else:
                                        merge_switch = False
                                        print(table_row.id, table_row._name, curr_repStart, curr_repEnd,table_row.repLeft,'no_merge', sep ='\t', file=open(log_filepath, 'a'))
                                        no_merge_frag_index.append(table_row._name)
                            else:
                                next_repName = table_by_id.iloc[row_num+1].repName
                                if next_repName != subfamily:
                                    merge_switch = False
                                    print(table_row.id, table_row._name, curr_repStart, curr_repEnd,table_row.repLeft,'end', sep ='\t', file=open(log_filepath, 'a'))
                                    merge_set.append(table_row._name)
                                    merge_frag_index.append(merge_set)

                                else:
                                    if curr_repLeft > common_length * 0.05:
                                        next_repStart = table_by_id.iloc[row_num+1].repStart
                                        next_repEnd = table_by_id.iloc[row_num+1].repEnd
                                        next_repLength = next_repEnd - next_repStart + 1
                                        overlap_covr = round((curr_repEnd - next_repStart+1)/genoLength*100,2)
                                        if (overlap_covr < overlap_cov_threshold) :
                                            #& (next_repLength <= common_length * 0.95)
                                            merge_switch = True
                                            print(table_row.id, table_row._name, curr_repStart, curr_repEnd,table_row.repLeft,'cont', overlap_covr, sep ='\t', file=open(log_filepath, 'a'))
                                            merge_set.append(table_row._name)        
                                        else:
                                            merge_switch = False
                                            print(table_row.id, table_row._name, curr_repStart, curr_repEnd,table_row.repLeft,'end', sep ='\t', file=open(log_filepath, 'a'))
                                            merge_set.append(table_row._name)
                                            merge_frag_index.append(merge_set)
                                    else:
                                        merge_switch = False
                                        print(table_row.id, table_row._name, curr_repStart, curr_repEnd,table_row.repLeft,'end', sep ='\t', file=open(log_filepath, 'a'))
                                        merge_set.append(table_row._name)
                                        merge_frag_index.append(merge_set)
                        else:
                            # last row
                            if not merge_switch:
                                merge_switch = False
                                print(table_row.id, table_row._name, curr_repStart, curr_repEnd,table_row.repLeft,'no_merge', sep ='\t', file=open(log_filepath, 'a'))
                                no_merge_frag_index.append(table_row._name)
                            else:
                                print(table_row.id, table_row._name, curr_repStart, curr_repEnd,table_row.repLeft,'end', sep ='\t', file=open(log_filepath, 'a'))
                                merge_set.append(table_row._name)
                                merge_frag_index.append(merge_set)
                    else:
                        # dealing with full length fragement
                        if not merge_switch:
                            merge_switch = False
                            print(table_row.id, table_row._name,curr_repStart, curr_repEnd,table_row.repLeft, 'no_merge' ,'complete', sep ='\t', file=open(log_filepath, 'a'))
                            no_mergr_complete_frag_index.append(table_row._name)
                        else:
                            merge_switch = False
                            print(table_row.id, table_row._name,curr_repStart, curr_repEnd,table_row.repLeft, 'end' ,'complete', sep ='\t', file=open(log_filepath, 'a'))
                            merge_set.append(table_row._name)
                            merge_frag_index.append(merge_set)        
                else:
                    print(table_row.id, table_row._name, table_row.repName, 'skip', sep ='\t', file=open(log_filepath, 'a'))
        total_frags_merged = 0
        for merge_set in merge_frag_index:
            total_frags_merged += len(merge_set)
        # final check
        # singlet = repeat entries with 1 fragment, no need to fix
        # all_complete_frag_index = repeat entries with multiple fragments, all full length
        # no_merge_frag_index = partial length fragments from repeat entries multiple fragments, no merge
        # no_mergr_complete_frag_index = full length fragments from repeat entries multiple fragments, no merge
        # total_frags_merged = merged fragments
        if subfam_table.shape[0] == (len(singlet_index)+len(all_complete_frag_index)+len(no_merge_frag_index)+len(no_mergr_complete_frag_index)+total_frags_merged):
            logging.info('checksum '+subfamily+': done')
        else:
            logging.info('checksum '+subfamily+': something is wrong, missing entry')
            next
        finalised_index_list = []
        internal_id_list = []
        internal_id = 0
        for idx in singlet_index:
            internal_id_list.append(internal_id)
            finalised_index_list.append(idx)
            internal_id += 1
        for idx in all_complete_frag_index:
            internal_id_list.append(internal_id)
            finalised_index_list.append(idx)
            internal_id += 1
        for idx in no_merge_frag_index:
            internal_id_list.append(internal_id)
            finalised_index_list.append(idx)
            internal_id += 1
        for idx in no_mergr_complete_frag_index:
            internal_id_list.append(internal_id)
            finalised_index_list.append(idx)
            internal_id += 1
        for merge_set in merge_frag_index:
            for idx in merge_set:
                internal_id_list.append(internal_id)
                finalised_index_list.append(idx)
            internal_id += 1
        
        dict_prep = {'rmsk_index': finalised_index_list,'internal_id': internal_id_list,}
        internal_id_tbl=pd.DataFrame(dict_prep)
        internal_id_tbl.to_csv(output_filepath, sep='\t')
        logging.info('done, saving the table at: '+output_filepath)
        logging.debug("--- %s seconds ---\n" % (time.time() - start_time))
#%%
overlap_cov_threshold = 100 
output_folder = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age'
assign_te_internal_id('L1MA10')
# %%
def main():
    subfamilies = config.subfamily_list
    global output_folder
    output_folder = config.internal_id_folder
    if os.path.isdir(output_folder) == False:
        os.mkdir(output_folder) 
    global overlap_cov_threshold
    overlap_cov_threshold = 100 # 100 = greedy, 0 = strict 
    with ProcessPoolExecutor(max_workers=40) as executor:
        results = executor.map(assign_te_internal_id,subfamilies)
    #for subfamily in subfamilies:
    #    assign_te_internal_id(subfamily)
#%%
if __name__ == '__main__':
    main()
