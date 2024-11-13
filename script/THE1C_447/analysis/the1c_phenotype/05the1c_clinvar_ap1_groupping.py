#%%
import sys
from turtle import width
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
from ma_mapper import mapper
from ma_mapper import custom_cmap
import numpy as np
#%%
subfamily = 'THE1C'
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
age_table = f'{config.te_age_folder}/{subfamily}.txt'
#%%
age_df = pd.read_csv(age_table, sep='\t')
internal_id_tbl = f'{config.internal_id_folder}/{subfamily}.txt'
internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
internal_id_sort = internal_id_df.sort_values('rmsk_index')
te_age_internal_id=internal_id_sort.merge(age_df, on='internal_id', how='left')
#%%
age_default_id = pd.DataFrame()
age_default_id['internal_id'] = subfamily + '_' + te_age_internal_id.index.astype(str)
age_default_id['te_age'] = te_age_internal_id['te_age']
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered, age_table=age_default_id)
#%%
coord_file=mapper.extract_metadata_from_alignment(alignment_file)
#%%
# %%
#field	example	description
#chrom	chr1	Chromosome (or contig, scaffold, etc.)
#chromStart	166070274	Start position in chromosome
#chromEnd	166070275	End position in chromosome
#name	name239343	Name of item
#score	3	P: 5, LP: 4, VUS: 3, LB: 2, B: 1, OTH: 0
#strand	.	+ or -
#thickStart	0	Start of where display should be thick (start codon)
#thickEnd	0	End of where display should be thick (stop codon)
#reserved	0,0,128	Used as itemRgb as of 2004-11-22
#lollySize	5	Size of lollipop
#changes	G>A(1)	changes
#variantIds	2344018	variantIds
#subIds	SCV003678877.2	subIds
#_mouseOver	chr1:166070275-166070275
#Variants (submissions):G>A(1)	mouseOver
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/phenotype_tracks/clinvarSubLolly.bed'
clinvar_df=pd.read_csv(bed_file, sep='\t', header = None)
clinvar_df.columns = ['chrom','start','end','name','score','strand','thickStart','thickEnd','reserved','lollySize','change','variantIds','subIds','_mouseOver']
#%%
clinvar_meat = clinvar_df[['chrom','start','end','name','score','strand']]
#%%
clinvar=mapper.map_and_overlay(alignment_file, coord_file, clinvar_meat, data_format='bed', custom_id=True, strand_overlap=False)

# %%
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/AP-1(bZIP).bed'
ap1=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed', custom_id=True, strand_overlap=True)
#%%
mean_clinvar=mapper.normalise(alignment=alignment_filtered, mapped_data=clinvar)
mean_ap1=mapper.normalise(alignment=alignment_filtered, mapped_data=ap1)
# %%
#find peaks
import scipy
peaks, _ = scipy.signal.find_peaks(mean_ap1, width = 6)
#highest_peak_index = peaks[np.argmax(mean_znf808[peaks])]
binding_indices = np.unique(np.where(ap1[:, peaks] != 0)[0])
nonbinding_indices=list(set(np.arange(ap1.shape[0])) - set(binding_indices))
#%%
# %%
clinvar_bind = clinvar[binding_indices]
clinvar_nonbind = clinvar[nonbinding_indices]
clinvar_sorted = np.vstack((clinvar_bind,clinvar_nonbind))
ap1_bind = ap1[binding_indices]
ap1_nonbind = ap1[nonbinding_indices]
ap1_sorted = np.vstack((ap1_bind,ap1_nonbind))
te_age_sorted=metadata_age.iloc[np.concatenate((binding_indices,nonbinding_indices))].te_age.fillna(0)
# %%
alignemnt_bind=alignment_filtered[binding_indices]
clinvar_bind[alignemnt_bind == 0] = np.nan
alignemnt_nonbind=alignment_filtered[nonbinding_indices]
clinvar_nonbind[alignemnt_nonbind == 0] = np.nan
stat_v_greater, p_value_greater = scipy.stats.mannwhitneyu(clinvar_bind,clinvar_nonbind, axis =0,nan_policy='omit', alternative = 'greater')
stat_v_less, p_value_less = scipy.stats.mannwhitneyu(clinvar_bind,clinvar_nonbind, axis =0,nan_policy='omit', alternative = 'less')

#%%
anno_label=metadata_age.fillna(0).te_age.sort_values().unique()
# %% 
from ma_mapper import plots
import importlib
importlib.reload(plots)
from ma_mapper import mapper
plots.plot_experimental(
    data = [clinvar_sorted,ap1_sorted,], 
    alignment=alignment_filtered, 
    heatmap_color=['clinvar','Greens'],
    heatmap_mode='spread_horizontal', 
    vlim = [[0,6],[0,0.0005],[0,0.0005]], 
    opacity = 0.5, 
    annotation=True, 
    anno_col = ['Blues'], 
    annotation_data=[te_age_sorted],
    anno_cbar_label=[anno_label],
    anno_title=['age'],
    anno_cbar_title=['MYA'], 
    aggregated=True, 
    aggregated_data=[mean_clinvar,-np.log10(p_value_greater),np.log10(p_value_less)], 
    agg_colset=['grey','blue','red'],
    agg_ylim=[[None,None]],
    agg_titles=['mean_phyloP','-log10P p>a','log10P p<a'], 
    #colorbar=True,
    #colorbar_steps = [0.1,0.0001],  
    )
#%%