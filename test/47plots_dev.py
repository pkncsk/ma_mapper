#%%
from configparser import Interpolation
import sys
from tabnanny import verbose
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import mapper
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import plots
from ma_mapper import custom_cmap
#%%
subfamily = ['THE1C']
alignment_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/alignment/'+subfamily[0]+'.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file)
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered, reference_table='/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/combined_age_div/combined_te_age_div.txt')
# %%
coord_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/coord_internal_id/'+subfamily[0]+'.txt'
#%%
bam_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/KZFP-bam_hg38/znf267.sorted.bam'
bam_forward=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_forward')
bam_reverse=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_reverse')
bam_min=mapper.map_and_overlay(alignment_file, coord_file, bam_file, data_format='read_min')
#%%
bigwig_file = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_bigwig/241-mammalian-2020v2.bigWig'
phylop=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig')
#%%
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/zoonomia447/hg38.phyloP447way.bw'
phylop_447=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig')
bigwig_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/zoonomia447/hg38.phyloP447wayLRT.bw'
phylop_LRT=mapper.map_and_overlay(alignment_file, coord_file, bigwig_file, data_format='bigwig')
#%%
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/AP-1(bZIP).bed'
ap1=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed')
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/NFkB-p65.bed'
nfkb=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed')
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/BMAL1.bed'
bmal=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed')
bed_file = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/_housekeeping/data/tfbs_homer/CLOCK.bed'
clock=mapper.map_and_overlay(alignment_file, coord_file, bed_file, data_format='bed')
#%%
#%% all overlay plot
import importlib
importlib.reload(plots)
importlib.reload(mapper)
plots.all_overlay_plot(alignment=alignment_filtered, alignment_col='nulc_white',data=[phylop_447, ap1],h_cmap=[custom_cmap.vlag_r_mpl, 'Greens'],vlim=[[-0.5,0.5],[0,1]], aggregated=True, a_colset = ['grey','green'], colorbar=True, cbar_steps=[0.1,0.1], heatmap_annot=metadata_age.te_age, opacity = 0.5, ag_ylim = [-0.5,0.5],logos=True)
#%% barchart only
plots.all_overlay_plot(alignment=alignment_filtered, alignment_col='nulc_white',data=[phylop_447, ap1, nfkb,bmal,clock], aggregated=True, a_colset = ['grey','green','red','blue','orange'], opacity = 0.5, heatmap=False, figsize=[20,10])
#%% barchart only
plots.all_overlay_plot(alignment=alignment_filtered, alignment_col='nulc_white',data=[phylop_447, bam_forward,bam_reverse], aggregated=True, a_colset = ['grey','red','blue','orange'], opacity = 0.5, heatmap=False, figsize=[20,10])
# %%
importlib.reload(plots)
importlib.reload(mapper)
plots.heatmap_overlay_agg_spread_plot(alignment=alignment_filtered, alignment_col='nulc_white',data=[phylop_447, ap1, nfkb,bmal,clock],h_cmap=[custom_cmap.vlag_r_mpl, 'Greens','Reds','Blues','Oranges'],vlim=[[-0.5,0.5],[0,1],[0,1],[0,1],[0,1]], aggregated=True, a_colset = ['grey','green','red','blue','orange'], colorbar=True, cbar_steps=[0.1,0.1,0.1,0.1,0.1], heatmap_annot=metadata_age.te_age, opacity = 0.5, ag_ylim = [-0.5,0.5],agg_titles=['phyloP','AP1','nfKb','BMAL','CLOCK'])
# %%
importlib.reload(plots)
importlib.reload(mapper)
plots.heatmap_spread_plot(alignment=alignment_filtered, alignment_col='dna',base_count=True,data=[phylop_447, ap1, nfkb,bmal,clock],h_cmap=[custom_cmap.vlag_r_mpl, 'Greens','Reds','Blues','Oranges'],vlim=[[-0.5,0.5],[0,1],[0,1],[0,1],[0,1]], aggregated=True, a_colset = ['grey','green','red','blue','orange'], colorbar=True, cbar_steps=[0.1,0.1,0.1,0.1,0.1], heatmap_annot=[metadata_age.te_age], opacity = 0.5, ag_ylim = [-0.5,0.5], spread_mode = 'vertical')
# %%
importlib.reload(plots)
importlib.reload(mapper)
plots.heatmap_spread_plot(alignment=alignment_filtered, alignment_col='dna',base_count=True,data=[phylop_447, ap1, nfkb,bmal,clock],h_cmap=[custom_cmap.vlag_r_mpl, 'Greens','Reds','Blues','Oranges'],vlim=[[-0.5,0.5],[0,1],[0,1],[0,1],[0,1]], aggregated=True, a_colset = ['grey','green','red','blue','orange'], colorbar=True, cbar_steps=[0.1,0.1,0.1,0.1,0.1], heatmap_annot=metadata_age.te_age, opacity = 0.5, ag_ylim = [-0.5,0.5])
# %%
alignment_filtered_complete, metadata_filtered= mapper.parse_and_filter(alignment_file, col_threshold = 0.5, col_content_threshold = 0.1, row_threshold = 0.5,)
from ma_mapper import sequence_alignment
importlib.reload(sequence_alignment)
sequence_alignment.parsed_array_to_sequence(alignment_filtered_complete,'THE1C','/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age/the1c_filtered_alignement.fasta')
# %%
