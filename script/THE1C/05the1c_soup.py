#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import extract_bigwig
from ma_mapper import extract_bam
from ma_mapper import mapper
#%%
subfamily = ['THE1C']
alignment_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/alignment/'+subfamily[0]+'.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered)
# %%
coord_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/coord_internal_id/'+subfamily[0]+'.txt'
#%%
bam_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/KZFP-bam_hg38/znf267.sorted.bam'
bam_forward=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='bam_forward')
bam_reverse=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='bam_reverse')
bam_min=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='bam_min')
#%%
bigwig_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_bigwig/241-mammalian-2020v2.bigWig'
phylop=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig')
#%%
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/AP-1(bZIP).bed'
ap1=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed')
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/NFkB-p65.bed'
nfkb=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed')
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/BMAL1.bed'
bmal=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed')
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/CLOCK.bed'
clock=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed')

# %%
from ma_mapper import plot
plot.overlay_plot([alignment_filtered,phylop,bam_forward, bam_reverse,ap1,nfkb,bmal,clock], metadata= metadata_age, nucleotide_color='white',data_vmin=[-0.5,0,0,0,0,0,0],data_vmax = [0.5,0.1,0.1,1,1,1,1],data_cmap=['RdBu','Purples','Oranges','Blues','Reds','Greens','copper'], plot_title='ChIPexo_znf276_tfbs', show_data_legend=True,data_label=['phyloP','ChIP_forward','ChIP_reverse','AP1','nfkB-p65','BMAL','CLOCK'])
# %%
import matplotlib.pyplot as plt
import numpy as np
normaliser = np.count_nonzero(alignment_filtered, axis=0)
fig, ax = plt.subplots(figsize=(10,3))
ax.fill_between(range(phylop.shape[1]), np.nansum(phylop, axis=0)/normaliser, color = 'grey')
ax.fill_between(range(bam_forward.shape[1]), y1= np.nansum(bam_forward, axis=0)/normaliser, color = 'purple',alpha = 0.5)
ax.fill_between(range(bam_reverse.shape[1]), np.nansum(bam_reverse, axis=0)/normaliser, color = 'orange',alpha = 0.5)
ax.fill_between(range(ap1.shape[1]), np.nansum(ap1, axis=0)/normaliser, color = 'blue',alpha = 0.5)
ax.fill_between(range(nfkb.shape[1]), np.nansum(nfkb, axis=0)/normaliser, color = 'red',alpha = 0.5)
ax.fill_between(range(bmal.shape[1]), np.nansum(bmal, axis=0)/normaliser, color = 'green',alpha = 0.5)
ax.fill_between(range(clock.shape[1]), np.nansum(clock, axis=0)/normaliser, color = 'yellow',alpha = 0.5)
ax.margins(x=0, y=0)
ax.set_ylim(-0.5,0.5)
# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
#peaks, peaks_info = find_peaks(np.nansum(bam_min, axis=0)/normaliser)
bam_min_post_process=np.minimum(np.nansum(bam_forward, axis=0)/normaliser,np.nansum(bam_reverse, axis=0)/normaliser)
global_max_idx=np.argmax(bam_min_post_process)
# %%
overlap_idx = []
nonoverlap_idx = []
for idx in range(len(bam_min)):
    if bam_min[idx][global_max_idx]> 0:
        overlap_idx.append(idx)
    else:
        nonoverlap_idx.append(idx)
overlap_phylop=phylop[overlap_idx,:]
nonoverlap_phylop=phylop[nonoverlap_idx,:]
#%%
normaliser_ov = np.count_nonzero(overlap_phylop, axis=0)
normaliser_nv = np.count_nonzero(nonoverlap_phylop, axis=0)
fig, ax = plt.subplots(figsize=(10,3))
ax.margins(x=0, y=0)
ax.plot(global_max_idx, 0, 'x')
ax.fill_between(range(overlap_phylop.shape[1]), np.nansum(overlap_phylop, axis=0)/normaliser_ov, color = 'blue',alpha=0.5)
ax.fill_between(range(nonoverlap_phylop.shape[1]), np.nansum(nonoverlap_phylop, axis=0)/normaliser_nv, color = ['red'],alpha = 0.5)
#%%
from scipy import stats
stat_v, p_value=stats.ttest_ind(overlap_phylop,nonoverlap_phylop, axis =0,nan_policy='omit')
fig, ax = plt.subplots(figsize=(10,3))
ax.plot(global_max_idx, 0, 'x')
ax.fill_between(range(len(p_value)), -np.log10(p_value), color='grey')
ax.fill_between(range(bam_forward.shape[1]), y1= bam_min_post_process, color = ['purple'],alpha = 0.5)
ax.hlines(y=2, xmin=0,xmax=len(p_value))
ax.margins(x=0, y=0)
# %%
from scipy import stats
stat_v, p_value=stats.ttest_ind(overlap_phylop,nonoverlap_phylop, axis =0,nan_policy='omit')
fig, ax = plt.subplots(figsize=(10,3))
ax.plot(global_max_idx, 0, 'x')
ax.fill_between(range(len(p_value)), -np.log10(p_value), color='grey')
ax.fill_between(range(bam_forward.shape[1]), y1= bam_min_post_process, color = ['purple'],alpha = 0.5)
ax.fill_between(range(ap1.shape[1]), np.nansum(ap1, axis=0)/normaliser, color = 'blue',alpha = 0.5)
ax.fill_between(range(nfkb.shape[1]), np.nansum(nfkb, axis=0)/normaliser, color = 'red',alpha = 0.5)
ax.fill_between(range(bmal.shape[1]), np.nansum(bmal, axis=0)/normaliser, color = 'green',alpha = 0.5)
ax.fill_between(range(clock.shape[1]), np.nansum(clock, axis=0)/normaliser, color = 'yellow',alpha = 0.5)
ax.hlines(y=2, xmin=0,xmax=len(p_value))
ax.margins(x=0, y=0)
# %%
from scipy import stats
stat_v, p_value=stats.mannwhitneyu(overlap_phylop,nonoverlap_phylop, axis =0,nan_policy='omit')
fig, ax = plt.subplots(figsize=(10,3))
ax.plot(global_max_idx, 0, 'x')
ax.fill_between(range(len(p_value)), -np.log10(p_value), color='grey')
ax.fill_between(range(bam_forward.shape[1]), y1= bam_min_post_process, color = ['purple'],alpha = 0.5)
ax.fill_between(range(ap1.shape[1]), np.nansum(ap1, axis=0)/normaliser, color = 'blue',alpha = 0.5)
ax.fill_between(range(nfkb.shape[1]), np.nansum(nfkb, axis=0)/normaliser, color = 'red',alpha = 0.5)
ax.fill_between(range(bmal.shape[1]), np.nansum(bmal, axis=0)/normaliser, color = 'green',alpha = 0.5)
ax.fill_between(range(clock.shape[1]), np.nansum(clock, axis=0)/normaliser, color = 'yellow',alpha = 0.5)
ax.hlines(y=2, xmin=0,xmax=len(p_value))
ax.margins(x=0, y=0)
# %%
from ma_mapper import extract_bed
matadata_bed=extract_bed.metadata_to_bed(coord_file, output_dir='/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/THE1C/the1c.bed',save_to_file=True, export_bedtool=True)
# %%
the1c_histone_enrichment=pd.read_csv('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/script/THE1C/the1c_enrichment.txt',sep='\t',header=None,names=['ID','Antigen class',	'Antigen',	'Cell class',	'Cell',	'Num of peaks',	'Overlaps / the1c',	'Overlaps / Control',	'Log P-val',	'Log Q-val',	'Fold Enrichment'])
#%%
the1c_histone_enrichment[['overlaps_exp','the1c']] =the1c_histone_enrichment['Overlaps / the1c'].str.split('/', expand=True)
# %%
the1c_activated=the1c_histone_enrichment[(the1c_histone_enrichment.Antigen=='H3K27ac') & (the1c_histone_enrichment['Cell class']=='Liver')]
# %%
the1c_activated=the1c_histone_enrichment[the1c_histone_enrichment.Antigen=='H3K27ac'].sort_values(['Log P-val','overlaps_exp'])
# %%
the1c_activated.sort_values(['overlaps_exp'])
# %%
bigwig_file='/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/chip_atlas/SRX5944575.bw'
cd4_h3k27ac=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig')
# %%
from ma_mapper import plot
plot.overlay_plot([alignment_filtered,cd4_h3k27ac], metadata= metadata_age, nucleotide_color='white',data_vmin=[0.00001,],data_vmax = [0.0001,],data_cmap=['Reds'], plot_title='the1c_h3k27ac', show_data_legend=True,data_label=['H3K27ac'])
# %%
plot.overlay_plot([alignment_filtered], metadata= metadata_age, nucleotide_color='white',plot_title='alignment')
# %%

#%%
blue_color_scale = [[0, 'rgba(0, 0, 0, 0.2)'],[1.0, 'rgba(8, 54, 116, 0.5)']]
five_color_scale = [[0, 'rgba(128, 128, 128, 0.2)'],[0.2, 'rgba(128, 128, 128, 0.2)'],[0.2, 'rgba(0, 128, 0, 0.2)'],[0.4, 'rgba(0, 128, 0, 0.2)'],[0.4, 'rgba(255, 255, 0, 0.2)'],[0.6, 'rgba(255, 255, 0, 0.2)'],[0.6, 'rgba(0, 0, 255, 0.2)'],[0.8, 'rgba(0, 0, 255, 0.2)'],[0.8, 'rgba(255, 0, 0, 0.2)'],[1.0, 'rgba(255, 0, 0, 0.2)']]
#%%
import plotly.graph_objects as go
subfamily = 'THE1C'
fig = go.Figure()
fig.add_trace(go.Heatmap(z= alignment_filtered,
                    colorscale=five_color_scale,
                    #text = text_array#colorbar=dict(ick0=0,dtick=1)
))
fig.add_trace(go.Heatmap(z= cd4_h3k27ac,
                    colorscale=blue_color_scale, zauto= False, zmax = 0.001
                    #colorbar=dict(tick0=0,dtick=1)
))
#fig.update_traces(text=text_array, texttemplate="%{text}")
fig.update_coloraxes(showscale=False)
fig.update_layout(xaxis_showgrid=False, xaxis_zeroline=False,violingap=0,violinmode='overlay',
    title=subfamily + ' Alignment',
    xaxis_title="Position (bp)",
    yaxis_title="Sequences",
    legend_title="base",
    showlegend=False,
    )
fig.update_layout(showlegend=False,)
fig.show()
# %%
# %%
overlap_idx = []
nonoverlap_idx = []
for idx in range(len(bam_min)):
    if bam_min[idx][global_max_idx]> 0:
        overlap_idx.append(idx)
    else:
        nonoverlap_idx.append(idx)
overlap_cd4_h3k27ac=cd4_h3k27ac[overlap_idx,:]
nonoverlap_cd4_h3k27ac=cd4_h3k27ac[nonoverlap_idx,:]
overlap_alignment = alignment_filtered[overlap_idx,:]
nonoverlap_alignment = alignment_filtered[nonoverlap_idx,:]
#%%
normaliser_ov = np.count_nonzero(overlap_alignment, axis=0)
normaliser_nv = np.count_nonzero(nonoverlap_alignment, axis=0)
fig, ax = plt.subplots(figsize=(10,3))
ax.margins(x=0, y=0)
ax.plot(global_max_idx, 0, 'x')
ax.fill_between(range(overlap_cd4_h3k27ac.shape[1]), np.nansum(overlap_cd4_h3k27ac, axis=0)/normaliser_ov, color = 'blue',alpha=0.5)
ax.fill_between(range(nonoverlap_cd4_h3k27ac.shape[1]), np.nansum(nonoverlap_cd4_h3k27ac, axis=0)/normaliser_nv, color = ['red'],alpha = 0.5)
# %%
from scipy import stats
stat_v, p_value=stats.mannwhitneyu(overlap_cd4_h3k27ac,nonoverlap_cd4_h3k27ac, axis =0,nan_policy='omit')
fig, ax = plt.subplots(figsize=(10,3))
ax.plot(global_max_idx, 0, 'x')
ax.fill_between(range(len(p_value)), -np.log10(p_value), color='grey')
ax.fill_between(range(bam_forward.shape[1]), y1= bam_min_post_process, color = ['purple'],alpha = 0.5)
ax.fill_between(range(ap1.shape[1]), np.nansum(ap1, axis=0)/normaliser, color = 'blue',alpha = 0.5)
ax.fill_between(range(nfkb.shape[1]), np.nansum(nfkb, axis=0)/normaliser, color = 'red',alpha = 0.5)
ax.fill_between(range(bmal.shape[1]), np.nansum(bmal, axis=0)/normaliser, color = 'green',alpha = 0.5)
ax.fill_between(range(clock.shape[1]), np.nansum(clock, axis=0)/normaliser, color = 'yellow',alpha = 0.5)
ax.hlines(y=2, xmin=0,xmax=len(p_value))
ax.margins(x=0, y=0)
# %%
cd4_h3k27ac_zeros = np.nan_to_num(cd4_h3k27ac)
rows_with_zeros = np.where(~cd4_h3k27ac_zeros.any(axis=1))[0]
rows_with_nonzeros = np.where(cd4_h3k27ac_zeros.any(axis=1))[0]


# %%
overlap_cd4_h3k27ac=phylop[rows_with_nonzeros,:]
nonoverlap_cd4_h3k27ac=phylop[rows_with_zeros,:]
overlap_alignment = alignment_filtered[rows_with_nonzeros,:]
nonoverlap_alignment = alignment_filtered[rows_with_zeros,:]
#%%
normaliser_ov = np.count_nonzero(overlap_alignment, axis=0)
normaliser_nv = np.count_nonzero(nonoverlap_alignment, axis=0)
fig, ax = plt.subplots(figsize=(10,3))
ax.margins(x=0, y=0)
ax.plot(global_max_idx, 0, 'x')
ax.fill_between(range(overlap_cd4_h3k27ac.shape[1]), np.nansum(overlap_cd4_h3k27ac, axis=0)/normaliser_ov, color = 'blue',alpha=0.5)
ax.fill_between(range(nonoverlap_cd4_h3k27ac.shape[1]), np.nansum(nonoverlap_cd4_h3k27ac, axis=0)/normaliser_nv, color = ['red'],alpha = 0.5)
# %%
from scipy import stats
stat_v, p_value=stats.mannwhitneyu(overlap_cd4_h3k27ac,nonoverlap_cd4_h3k27ac, axis =0,nan_policy='omit')
fig, ax = plt.subplots(figsize=(10,3))
ax.plot(global_max_idx, 0, 'x')
ax.fill_between(range(len(p_value)), -np.log10(p_value), color='grey')
ax.fill_between(range(bam_forward.shape[1]), y1= bam_min_post_process, color = ['purple'],alpha = 0.5)
ax.fill_between(range(ap1.shape[1]), np.nansum(ap1, axis=0)/normaliser, color = 'blue',alpha = 0.5)
ax.fill_between(range(nfkb.shape[1]), np.nansum(nfkb, axis=0)/normaliser, color = 'red',alpha = 0.5)
ax.fill_between(range(bmal.shape[1]), np.nansum(bmal, axis=0)/normaliser, color = 'green',alpha = 0.5)
ax.fill_between(range(clock.shape[1]), np.nansum(clock, axis=0)/normaliser, color = 'yellow',alpha = 0.5)
ax.hlines(y=2, xmin=0,xmax=len(p_value))
ax.margins(x=0, y=0)
#%%