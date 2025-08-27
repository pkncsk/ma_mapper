#%%
from ma_mapper import mapper
from ma_mapper import plots
from ma_mapper import custom_cmap
import numpy as np
import pandas as pd
#%%

subfamily='THE1C'
alignment_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/alignment/THE1C.fasta.aligned'
genomewide_data_filepath  = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/UCSC_phyloP_track/hg38.phyloP447way.bw'
#%%

#%%
filtered_alignment_matrix, alignment_coordinate  = mapper.parse_and_filter(alignment_file=alignment_filepath,col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10,)
data_matrix=mapper.map_and_overlay(alignment_filepath, genomewide_data_filepath,data_format='bigwig', col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
#%%
bed_filepath = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/peaks/kzfp_peak_bed/hg38_kzfps_combined.bed'
kzfp_df=pd.read_csv(bed_filepath, sep='\t', header=None)
kzfp_df.columns=['chrom','start','end','name','score','strand']
bed_file=kzfp_df[kzfp_df['name'].str.contains('ZNF267')]
znf267=mapper.map_and_overlay(alignment_filepath, bed_file, data_format='bed', strand_overlap=False, col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
#%%
atac_df_bed = pd.read_csv('/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/phenotype_tracks/atac_coord_hg38.bed', sep='\t', header = None)
atac=mapper.map_and_overlay(alignment_filepath, atac_df_bed, data_format='bed', strand_overlap=False, col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
#%%
from matplotlib import cm
from matplotlib.colors import ListedColormap

brighter_cmap = cm.get_cmap('RdBu')(np.linspace(0.05, 0.95, 256))
from ma_mapper import plots
from ma_mapper import custom_cmap
plots.plot(
    data = [znf267,atac], 
    alignment=filtered_alignment_matrix,
    show_alignment=False, 
    heatmap_color=['Purples','Yellows'],
    heatmap_mode='overlay',
    heatmap_title=['ZNF267 and ATAC-seq binding signals on THE1C MSA'],
    heatmap_title_fs = 10, 
    heatmap_xlabel = 'position (bp)',
    heatmap_xlabel_fs = 10,
    heatmap_ylabel = 'sequences',
    heatmap_ylabel_fs = 10,
    vlim = [[0,7],[0,1]], 
    opacity = 0.8, 
    hm_transparency_mode= 'static',
    hm_interpolation='auto',
    colorbar=False,
    image_res=300
    )
#%%
from ma_mapper import plots
from ma_mapper import custom_cmap
plots.plot(
    data = [data_matrix], 
    alignment=filtered_alignment_matrix,
    show_alignment=False, 
    heatmap_color=['RdBu'],
    heatmap_mode='overlay',
    heatmap_title=['PhyloP from Zoonomia dataset on THE1C MSA'],
    heatmap_title_fs = 10, 
    heatmap_xlabel = 'position (bp)',
    heatmap_xlabel_fs = 10,
    heatmap_ylabel = 'sequences',
    heatmap_ylabel_fs = 10,
    vlim = [[-1,1]], 
    opacity = 1.0, 
    hm_transparency_mode= 'gradient',
    colorbar=True,
    colorbar_steps=[0.1],
    xlim=[100,150]
    )
# %%
from ma_mapper import plots
from ma_mapper import custom_cmap
plots.plot(
    alignment=filtered_alignment_matrix,
    show_alignment=True, 
    show_alignment_colbar=True,
    alignment_col='dna',
    heatmap_mode='overlay',
    heatmap_title=['THE1C MSA'],
    heatmap_title_fs = 10, 
    heatmap_xlabel = 'position (bp)',
    heatmap_xlabel_fs = 10,
    heatmap_ylabel = 'sequences',
    heatmap_ylabel_fs = 10,
    opacity = 1.0, 
    hm_transparency_mode= 'gradient',
    )
# %%
plots.plot(
    alignment=filtered_alignment_matrix,
    show_alignment=True, 
    show_alignment_colbar=True,
    alignment_col='dna',
    heatmap_mode='overlay',
    heatmap_title=['THE1C MSA at 200-250bp'],
    heatmap_title_fs = 10, 
    heatmap_xlabel = 'position (bp)',
    heatmap_xlabel_fs = 10,
    heatmap_ylabel = 'sequences',
    heatmap_ylabel_fs = 10,
    opacity = 1.0, 
    hm_transparency_mode= 'gradient',
    xlim=[200,250]
    )
# %%
filtered_alignment_matrix, alignment_coordinate  = mapper.parse_and_filter(alignment_file=alignment_filepath,col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10, extension_length=100,source_fasta='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/human_genome_fasta/hg38_fasta/hg38.fa')
data_matrix=mapper.map_and_overlay(alignment_filepath, genomewide_data_filepath,data_format='bigwig', col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10, extension_length=100,source_fasta='/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/data/resource/human_genome_fasta/hg38_fasta/hg38.fa')
#%%
mean_phylop=mapper.normalise(alignment_matrix=filtered_alignment_matrix, data_matrix=data_matrix)
mean_track=mean_phylop[100:-100]
consensus_length=len(mean_track)
pos_posit=len(mean_track[mean_track>0])
neg_posit=len(mean_track[mean_track<0])
#%%
from matplotlib import cm
from matplotlib.colors import ListedColormap
brighter_cmap = cm.get_cmap('RdBu')(np.linspace(0.1, 0.9, 256))
from ma_mapper import plots
from ma_mapper import custom_cmap
plots.plot(
    data = [data_matrix], 
    alignment=filtered_alignment_matrix,
    show_alignment=False, 
    heatmap_color=[ListedColormap(brighter_cmap)],
    heatmap_mode='overlay',
    heatmap_title=['Zoonomia phyloP of THE1C MSA with 100 bp flanks'],
    heatmap_title_fs = 10, 
    heatmap_xlabel = 'position (bp)',
    heatmap_xlabel_fs = 10,
    heatmap_ylabel = 'sequences',
    heatmap_ylabel_fs = 10,
    vlim = [[-1,1]], 
    colorbar=True,
    colorbar_steps=[0.1],
    agg_major_tick = 200,
    aggregated=True, 
    aggregated_data=[mean_phylop], 
    agg_colset=['grey',],
    agg_ylim=[[-1.2,1.2]],
    #agg_ylabel_right=['ap1 motif'], 
    agg_ylabel=['mean phyloP'],
    agg_ylabel_fs=10,
    agg_xlabel='position (bp)',
    hm_interpolation='auto',
    agg_plottext=[f'consensus length: {consensus_length}\n phyloP>0: {pos_posit}\n phyloP<0: {neg_posit}'],
    agg_plottext_fs=6,
    agg_plottext_pos=[0.99,0.65],
    save_to_file= 'phyloP_zoomed.png',
    )
# %%
filtered_alignment_matrix, alignment_coordinate  = mapper.parse_and_filter(alignment_file=alignment_filepath,col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10,)
phylop_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/UCSC_phyloP_track/hg38.phyloP447way.bw'
phylop_matrix=mapper.map_and_overlay(alignment_filepath, genomewide_data_filepath,data_format='bigwig', col_threshold = 0.10, col_content_threshold = 0.10, row_threshold = 0.10)
mean_phylop=mapper.normalise(alignment_matrix=filtered_alignment_matrix, data_matrix=phylop_matrix)
#%%
gnomad_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/allele_frequency_vcf/gnomad'
gnomad_matrix=mapper.map_and_overlay(alignment_filepath, gnomad_filepath,data_format='vcf', vcf_format='gnomad', pf_col_threshold = 0.10, pf_col_content_threshold = 0.10, pf_row_threshold = 0.10)
mean_gnomad=mapper.normalise(alignment_matrix=filtered_alignment_matrix, data_matrix=gnomad_matrix)
# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

# Reshape for sklearn
X = np.array(mean_phylop).reshape(-1, 1)
y = np.array(mean_gnomad)

# Fit regression model
model = LinearRegression()
model.fit(X, y)
y_pred = model.predict(X)

# R² and adjusted R²
r2 = r2_score(y, y_pred)
n = len(y)
p = 1  # one predictor (phyloP)
adj_r2 = 1 - (1 - r2) * (n - 1) / (n - p - 1)

# Pearson and Spearman
pearson_r, pearson_p = stats.pearsonr(mean_phylop, mean_gnomad)
spearman_rho, spearman_p = stats.spearmanr(mean_phylop, mean_gnomad)

# Plot
fig, ax = plt.subplots(figsize=(8, 8))
ax.scatter(X, y, alpha=0.5, color='black', marker='.')

# Plot regression line
ax.plot(X, y_pred, color='red', linewidth=2, alpha=0.7)

# Annotate
ax.text(
    0.99, 0.80,
    f"n = {n}\n"
    f"Pearson r = {pearson_r:.2f} (p = {pearson_p:.2e})\n"
    f"Spearman ρ = {spearman_rho:.2f} (p = {spearman_p:.2e})\n"
    f"y = {model.coef_[0]:.2e}x + {model.intercept_:.2e}\n"
    f"R² = {r2:.2f}, adj. R² = {adj_r2:.2f}",
    transform=ax.transAxes,
    fontsize=11,
    verticalalignment='top',
    horizontalalignment='right',
    bbox=dict(facecolor='white', alpha=0.85, edgecolor='gray')
)

# Axes styling
ax.axvline(x=0, color='grey', linewidth=1, alpha=0.5)
ax.set_ylim(0, 0.003)
ax.set_xlabel('Mean phyloP per base', fontsize=12)
ax.set_ylabel('Mean alternate allele frequency per base', fontsize=12)
ax.set_title('Mean phyloP vs. mean alternate allele frequency\nper base on THE1C MSA', fontsize=14)

plt.tight_layout()
plt.savefig('./compare.png', dpi=300, bbox_inches="tight")
plt.show()

# %%
