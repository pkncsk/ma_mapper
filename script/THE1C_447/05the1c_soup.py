#%%
import sys
import pandas as pd
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import custom_cmap
from ma_mapper import mapper
#%%
subfamily = 'THE1C'
alignment_file = f'{config.te_alignment_folder}/{subfamily}.fasta.aligned'
alignment_filtered, metadata_filtered= mapper.parse_and_filter(alignment_file,)
age_table = f'{config.te_age_folder}/{subfamily}.txt'
age_df = pd.read_csv(age_table, sep='\t')
internal_id_tbl = f'{config.internal_id_folder}/{subfamily}.txt'
internal_id_df = pd.read_csv(internal_id_tbl, sep='\t')
te_age_internal_id=internal_id_df.merge(age_df, on='internal_id')
age_default_id = pd.DataFrame()
age_default_id['internal_id'] = subfamily + '_' + te_age_internal_id.index.astype(str)
age_default_id['te_age'] = te_age_internal_id['te_age']
metadata_age = mapper.match_age_to_id_metadata(metadata_filtered, age_table=age_table)
# %%
coord_file=mapper.extract_metadata_from_alignment(alignment_file)
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
mean_phylop=mapper.normalise(alignment=alignment_filtered, mapped_data=phylop_447)
mean_ap1=mapper.normalise(alignment=alignment_filtered, mapped_data=ap1)
mean_nfkb=mapper.normalise(alignment=alignment_filtered, mapped_data=nfkb)
mean_bmal=mapper.normalise(alignment=alignment_filtered, mapped_data=bmal)
mean_clock=mapper.normalise(alignment=alignment_filtered, mapped_data=clock)
#%%
from ma_mapper import plots
import importlib
importlib.reload(plots)
plots.plot_experimental(alignment=alignment_filtered, alignment_col='nulc_white',data=[phylop_447, ap1, nfkb,bmal,clock],h_cmap=[custom_cmap.vlag_r_mpl, 'Greens','Reds','Blues','Oranges'],vlim=[[-0.5,0.5],[0,1],[0,1],[0,1],[0,1]], aggregated=True, a_colset = ['grey','green','red','blue','orange'],aggregated_data=[mean_phylop, mean_ap1, mean_nfkb, mean_bmal, mean_clock], colorbar=True, cbar_steps=[0.1,0.1,0.1,0.1,0.1], annotation=True,annotation_data=[metadata_age.te_age.fillna(0)],anno_col=['Greens'], opacity = 0.5, agg_ylim = [None],agg_titles=['phyloP','AP1','nfKb','BMAL','CLOCK'])
#%%