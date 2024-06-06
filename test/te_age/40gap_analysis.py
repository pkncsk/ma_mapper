#%%
from functools import total_ordering
from itertools import accumulate
from mailcap import getcaps
import pandas as pd
import numpy as np
import sys
import os
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/stable')
import config_hg38_repeatlib as config
from concurrent.futures import ProcessPoolExecutor
#%%
combined_te_gap=pd.read_csv(f'/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/gap/combined_te_gap.txt', sep='\t')
# %%
# %%
histo_gap=combined_te_gap.groupby('gap').count()
# %%
3284146/4624576 standalone/total (71% of all TE entries from repeatmasker)
1340430 fragments (29%) were joined by algorithm 
547387 joined fragments from 793043 gaps 
gap breakdown
gap/overlap         count   percent accumulated perc.
overlaps            33769   4.26    4.26
0                   7819    0.98    5.24     
0<length<10         19609   2.47    7.71 
10<length<=100      109816  13.85   21.56
100<length<=500     431033  54.35   75.91 
500<length<=1000    90670   11.43   87.34
1000<length<=5000   86781   10.94   98.28
5000<length<=10000  12322   1.55    99.83
length>10000        1224    0.15    99.97
total gap           793043  100
#%%
id_count=combined_te_gap.fillna(0).groupby('internal_id').count()
#%%
combined_te_gap.gap.max()
combined_te_gap[combined_te_gap.gap==31384]
# %%
subfamily = 'L1MA10'
subfamily_filename = subfamily.replace('/','_')
filtered_table = config.filtered_table
internal_id_folder = config.internal_id_folder
internal_id_tbl = f'{internal_id_folder}/{subfamily_filename}.txt'
internal_id_df = pd.read_csv(internal_id_tbl, sep='\t', index_col= 0)
subfam_df = filtered_table[filtered_table.repName == subfamily]
subfam_w_internal_id=pd.merge(subfam_df, internal_id_df, left_index = True, right_on = 'rmsk_index')
# %%
