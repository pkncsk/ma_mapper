#%% LOAD PACKAGE
from ma_mapper import mapper
from ma_mapper import plots
from ma_mapper import custom_cmap
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
def fdr_control_with_nans(p_values, method='bh'):
    p_values = np.asarray(p_values)
    valid_mask = ~np.isnan(p_values)
    valid_p_values = p_values[valid_mask]
    adj_p_values_valid = false_discovery_control(valid_p_values, method=method)
    adj_p_values = np.full(p_values.shape, np.nan)
    adj_p_values[valid_mask] = adj_p_values_valid
    return adj_p_values
#%% INITIAL PARAMETER
subfamily='THE1C'
alignment_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/alignment/THE1C.fasta.aligned'
genomewide_data_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/annotation/homer_known_motif_hg38/BMAL1(bHLH).bed'
phyloP_data_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/UCSC_phyloP_track/hg38.phyloP447way.bw'
te_age_folder = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime/age_lenient/'
age_table = f'{te_age_folder}/{subfamily}.txt'
internal_id_folder = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/output/teatime/internal_id'
internal_id_tbl = f'{internal_id_folder}/{subfamily}.internal_id.txt'
#%% SIMPLE WORKFLOW
filtered_alignment_matrix, alignment_coordinate  = mapper.parse_and_filter(alignment_file=alignment_filepath,col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
motif_matrix=mapper.map_and_overlay(alignment_filepath, genomewide_data_filepath,data_format='bed',strand_overlap =True, pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
phyloP_matrix=mapper.map_and_overlay(alignment_filepath, phyloP_data_filepath,data_format='bigwig', pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
age_df = pd.read_csv(age_table, sep='\t')
internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
internal_id_sort = internal_id_df.sort_values('rmsk_index')
te_age_internal_id=internal_id_sort.merge(age_df, on='internal_id', how='left')
#%%
age_default_id = pd.DataFrame()
age_default_id['internal_id'] = subfamily + '_' + te_age_internal_id.index.astype(str)
age_default_id['te_age'] = te_age_internal_id['te_age']
#%%
coord_age = mapper.match_age_to_id_coordinate(alignment_coordinate, age_table=age_default_id)
###############IMPORTANT####################
#filter NA
coord_age=coord_age[~coord_age['te_age'].isna()]
noNA_indices=coord_age.index
coord_age=coord_age.reset_index()
filtered_alignment=filtered_alignment_matrix[noNA_indices]

coord_age['len'] = coord_age.end.astype(int) - coord_age.start.astype(int)
age_subgroups = np.unique(coord_age['te_age'].sort_values())
age_subgroup = {subgroup: num for num, subgroup in enumerate(age_subgroups)}
age_anno=coord_age['te_age'].map(age_subgroup)
# %%
phyloP_matrix=phyloP_matrix[noNA_indices]
motif_matrix=motif_matrix[noNA_indices]
mean_phylop=mapper.normalise(alignment_matrix=filtered_alignment, data_matrix=phyloP_matrix)
cov_motif_matrix=mapper.normalise(alignment_matrix=filtered_alignment, data_matrix=motif_matrix, method='perc_coverage')
#%%
from ma_mapper import plots
from ma_mapper import mapper
plots.plot(
    data = [motif_matrix,], 
    alignment=filtered_alignment,
    show_alignment=True, 
    heatmap_color=['Blues'],
    heatmap_mode='overlay',
    heatmap_title=['NFkB motifs on THE1C MSA'],
    heatmap_title_fs = 10, 
    anno_ylabel = 'sequences',
    anno_ylabel_fs=10,
    vlim = [[0,1]], 
    opacity = 0.9, 
    hm_transparency_mode= 'gradient',
    aggregated=True, 
    aggregated_data=[cov_motif_matrix], 
    agg_colset=['blue',],
    agg_ylim=[[0,50]],
    #agg_ylabel_right=['ap1 motif'], 
    agg_ylabel=['perc_coverage'],
    agg_ylabel_fs=7,
    agg_xlabel='position (bp)',
    annotation=False, 
    anno_col = ['Blues'], 
    annotation_data=[age_anno.values],
    anno_cbar=True,
    anno_cbar_label=[age_subgroups],
    anno_cbar_title=['TEA-TIME'], 
    colorbar=False,

    )
#%%
# %%
import scipy
peaks, _ = scipy.signal.find_peaks(cov_motif_matrix, width = 6)
#%%
#right,left
#highest_peak_index = peaks[np.argmax(mean_znf808[peaks])]
binding_indices = np.unique(np.where(motif_matrix[:, peaks] != 0)[0])
nonbinding_indices=list(set(np.arange(motif_matrix.shape[0])) - set(binding_indices))
#%%
binding_indices_right = np.unique(np.where(motif_matrix[:, 180 ] != 0)[0])
nonbinding_indices_right=list(set(np.arange(motif_matrix.shape[0])) - set(binding_indices_right))
binding_indices_left = np.unique(np.where(motif_matrix[:, 115] != 0)[0])
nonbinding_indices_left=list(set(np.arange(motif_matrix.shape[0])) - set(binding_indices_left))
#%%
coord_age['right_group'] = 'No Group'
coord_age.loc[coord_age.index.isin(binding_indices_right), 'right_group'] = 'A'
coord_age.loc[coord_age.index.isin(nonbinding_indices_right), 'right_group'] = 'B'
coord_age['left_group'] = 'No Group'
coord_age.loc[coord_age.index.isin(binding_indices_left), 'left_group'] = 'A'
coord_age.loc[coord_age.index.isin(nonbinding_indices_left), 'left_group'] = 'B'


#%%
coord_sorted=coord_age.sort_values(['right_group','left_group','te_age'])
sorted_indices=coord_sorted.index
subgroups = np.unique(coord_sorted['right_group'].astype(str))
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno_right=coord_sorted['right_group'].map(numerical_subgroup)
subgroups = np.unique(coord_sorted['left_group'].astype(str))
numerical_subgroup = {subgroup: num for num, subgroup in enumerate(subgroups)}
subgroup_anno_left=coord_sorted['left_group'].map(numerical_subgroup)
#%%
phylop_sorted=phyloP_matrix[sorted_indices]
motif_matrix_sorted=motif_matrix[sorted_indices]
alignment_sorted=filtered_alignment[sorted_indices]
age_anno_sorted=age_anno[sorted_indices]
#%%
binding_indices=coord_sorted[(coord_sorted['right_group']=='A')|(coord_sorted['left_group']=='A')].index
nonbinding_indices=coord_sorted[(coord_sorted['right_group']=='B')&(coord_sorted['left_group']=='B')].index
only_right = coord_sorted[(coord_sorted['right_group']=='A')&(coord_sorted['left_group']=='B')].index
only_left = coord_sorted[(coord_sorted['right_group']=='B')&(coord_sorted['left_group']=='A')].index
intersect_right_left = coord_sorted[(coord_sorted['right_group']=='A')&(coord_sorted['left_group']=='A')].index
#%%
binding_indices=coord_age[(coord_age['right_group']=='A')|(coord_age['left_group']=='A')].index
nonbinding_indices=coord_age[(coord_age['right_group']=='B')&(coord_age['left_group']=='B')].index
phylop_bind = phyloP_matrix[binding_indices]
phylop_nonbind = phyloP_matrix[nonbinding_indices]
alignemnt_bind=filtered_alignment[binding_indices]
phylop_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=filtered_alignment[nonbinding_indices]
phylop_nonbind[alignemnt_nonbind == 0] = np.nan
alignment_sorted = np.vstack((alignemnt_bind,alignemnt_nonbind))
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
p_value_greater_adjusted=fdr_control_with_nans(p_value_greater)
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(phylop_bind,phylop_nonbind, axis =0,nan_policy='omit', alternative = 'less')
p_value_less_adjusted=fdr_control_with_nans(p_value_less)
# %%
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot(
    heatmap=False,
    show_alignment=False,
    alignment=alignment_sorted,
    #logos=True, 
    aggregated=True, 
    aggregated_data=[-np.log10(p_value_greater_adjusted),np.log10(p_value_less_adjusted)], 
    agg_colset=['blue','red'],
    agg_ylim=[[0,10],[-10,0]],
    agg_yhighlight=[[109,122],[175,190],],
    agg_yhighlight_col= ['blue','blue'],
    agg_yhighlight_alpha=[0.2,0.2,0.2,0.2],
    agg_xhighlight=[[3,3],[-3,-3]],
    agg_xhighlight_col= ['grey','grey'],
    agg_xhighlight_alpha=[0.2,0.2],
    agg_ylabel_right=['conserved\nmotif>no motif','accelerated\nmotif<no motif'], 
    agg_ylabel_right_fs=10,
    agg_ylabel_right_pos=[1.2,0.5],
    agg_ylabel=['-log10(adjusted p-value)',None],
    agg_ylabel_ypos=[0.05,None],
    agg_ylabel_xpos=0,
    agg_ylabel_fs=10,
    #agg_plot_title=['phyloP of THE1C, grouped by AP-1 motif',None],
    agg_plot_title_fs= 12, 
    agg_xlabel = 'position (bp)',
    agg_xlabel_fs=10,
    agg_plottext=[None,f'TE w/ motifs: {len(binding_indices)}\nTE w/o motifs: {len(nonbinding_indices)}'],
    agg_plottext_fs=10,
    agg_plottext_pos=[0.99,0.01],
    agg_major_tick=20,
    colorbar=False,
    agg_h=20,
    figsize= [80,40],
    #xlim=[304,313],
    #gg_major_tick=1,
    )
#%%
binding_alignment=filtered_alignment_matrix[binding_indices]
nonbinding_alignment = filtered_alignment_matrix[nonbinding_indices]
# %%
motif_range=binding_alignment
#%%
motif_range = binding_alignment[:, 25:75]
#%%
#motif_range = nonbinding_alignment[:, 45:55]
#%%
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
base_mapping = {0: 'N', 1: 'A', 2: 'C', 3: 'T',4:'G'}
name = 'entry'
# Convert each row to DNA sequence
dna_sequences = []
for idx, row in enumerate(motif_range):
    seq_name = f'{name}_{idx}'
    dna_sequence = ''.join(base_mapping[int(val)] for val in row)
    dna_sequences.append(SeqRecord(Seq(dna_sequence),seq_name , '', ''))

import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# --- Step 1: Slice motif region from filtered alignment matrix ---  # shape: (num_seqs, 21)

# --- Step 2: Convert alignment matrix to SeqRecord list ---
base_mapping = {0: 'N', 1: 'A', 2: 'C', 3: 'T', 4: 'G'}
name_prefix = 'entry'

all_records = []
for idx, row in enumerate(motif_range):
    seq_str = ''.join(base_mapping[int(val)] for val in row)
    rec = SeqRecord(Seq(seq_str), id=f'{name_prefix}_{idx}', description='')
    all_records.append(rec)

# --- Step 3: Filter out sequences with ambiguous bases or gaps ---
cleaned_records = [rec for rec in all_records if 'N' not in str(rec.seq).upper() and '-' not in str(rec.seq)]

# --- Step 4: Remove exact duplicates ---
unique_seq_dict = {}
for rec in cleaned_records:
    seq_str = str(rec.seq).upper()
    unique_seq_dict[seq_str] = rec  # Keep only one representative per sequence

non_redundant_records = list(unique_seq_dict.values())

# --- Step 5: Save to FASTA for MEME ---
SeqIO.write(non_redundant_records, "non_redundant_for_meme.fasta", "fasta")

import numpy as np
import pandas as pd
import logomaker
import matplotlib.pyplot as plt

# --- STEP 0: Input ---
# Example matrix (your actual input will be similar)
# matrix = np.load("your_matrix.npy")  # (N, L)
 # Simulated example
matrix=motif_range
# --- STEP 1: Filter out rows with '0' (ambiguous base N) ---
clean_matrix = matrix[~np.any(matrix == 0, axis=1)]
print(f"Matrix shape after removing Ns: {clean_matrix.shape}")

# --- STEP 2: Count base frequency at each position ---
# Base index map: 1=A, 2=C, 3=T, 4=G → offset to 0=A, 1=C, 2=G, 3=T
# We'll map to rows ['A', 'C', 'G', 'T']
def count_bases(mat):
    counts = np.zeros((4, mat.shape[1]), dtype=int)
    for i, base in enumerate(['A', 'C', 'T', 'G']):
        # base code: 1=A, 2=C, 3=T, 4=G → shift by 1
        counts[i] = np.sum(mat == (i + 1), axis=0)
    return pd.DataFrame(counts, index=['A', 'C', 'T', 'G'])

counts_df = count_bases(clean_matrix)
print("Raw base counts:")
print(counts_df)

# --- STEP 3: Convert to PWM (Position Probability Matrix) ---
ppm_df = logomaker.transform_matrix(counts_df.T, from_type='counts', to_type='probability').T

# --- STEP 4: Visualize the motif with Logomaker ---
plt.figure(figsize=(6, 2))
logomaker.Logo(ppm_df.T)
plt.title("Sequence Logo (PWM)")
plt.savefig('motiflogos.png', dpi=300)
plt.tight_layout()
plt.show()

# --- STEP 5: Write to MEME motif file ---
def write_meme_motif(ppm_df, motif_name, outfile, nsites=None, e_value=0.0):
    with open(outfile, 'w') as f:
        f.write("MEME version 4\n\n")
        f.write("ALPHABET= ACGT\n\n")
        f.write("strands: + -\n\n")
        f.write("Background letter frequencies\n")
        f.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")
        f.write(f"MOTIF {motif_name}\n")
        f.write(f"letter-probability matrix: alength= 4 w= {ppm_df.shape[1]} nsites= {nsites or clean_matrix.shape[0]} E= {e_value}\n")

        # Reorder to A, C, G, T as expected by MEME format
        for i in range(ppm_df.shape[1]):
            col = ppm_df.iloc[:, i]
            f.write(" ".join(f"{col[base]:.6f}" for base in ['A','C','G','T']) + "\n")

# --- STEP 6: Export ---
write_meme_motif(ppm_df, motif_name="TE_alignment_motif", outfile="custom_pwm.meme")
print("✅ MEME motif file saved as: custom_pwm.meme")

# %%
