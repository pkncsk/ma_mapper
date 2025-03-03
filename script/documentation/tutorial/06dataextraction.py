#%% load package
from ma_mapper import mapper

alignment_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/ma_mapper/hg38_main/alignment/THE1C.fasta.aligned'
BED_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/annotation/homer_known_motif_hg38/AP-1(bZIP).bed'
#extract TE coordinates from alignment fasta the header for the alignment should be in >NAME::chromosome:start-stop(strand) format
alignment_matrix, alignment_coordinate, filters  = mapper.parse_and_filter(alignment_file=alignment_filepath, preprocess_out=True)
#%% BED file 
# even though mapper module can handle general tasks, we can use dedicated modules for data extraction, with extra steps
# strating with BED file which was showned briefly in the basicworkflow section
# extract_bed is a wrapper for pybedtools
from ma_mapper import extract_bed
#then use alignment coordinate to extract genomewide data
output_matrix=extract_bed.bed_io(coordinate_table=alignment_coordinate, bed=BED_filepath)

#the default setting for BED file extraction is to match DNA strand, in the case where strand matching is not needed, it can be disabled
output_matrix=extract_bed.bed_io(coordinate_table=alignment_coordinate, bed=BED_filepath, strand_overlap=False)
#%% BIGWIG file 
# was also introduced earlier in basic workflow, as of now there is no options for this function
# extract_bigwig is a wrapper for pybigwig
from ma_mapper import extract_bigwig
BIGWIG_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/UCSC_phyloP_track/hg38.phyloP447way.bw'
extract_bigwig.bigwig_io(coordinate_table=alignment_coordinate, bigwig=BIGWIG_filepath)
#%% BAM/SAM file
# extract_bam is a wrapper for pysam
from ma_mapper import extract_bam
BAM_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/sequenceing_data_fastqsambam/znf267.sorted.bam'
# this module has different outputs which can be extracted, we can separate reverse, forward reads or have them all as summation, and have all those as smooted reads
read_forward_matrix=extract_bam.bam_io(coordinate_table=alignment_coordinate,bam_file=BAM_filepath, bam_format='read_forward')
read_reverse_matrix=extract_bam.bam_io(coordinate_table=alignment_coordinate,bam_file=BAM_filepath, bam_format='read_reverse')
read_total_matrix=extract_bam.bam_io(coordinate_table=alignment_coordinate,bam_file=BAM_filepath, bam_format='read_sum')
# it is also possible to control probe and smoothing length for smooted reads
normal_forward_matrix=extract_bam.bam_io(coordinate_table=alignment_coordinate,bam_file=BAM_filepath, bam_format='normal_forward', probe_length=100, smoothing_length=100)
#%% MAF file
# extract_maf is a wrapper for MafIO from biopython package
# it also modified get_spliced function to handle duplicate entries in zoonomia maf file, as for now it will collaspe duplicates according to IUPAC rules, add parallelization support to speed up the process
from ma_mapper import extract_maf
MAF_dir = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/multi_species_multiple_alignment_maf/cactus447'
# since maf is a collection of multiple seqeucne alignment, extracting the data into numerical matrix means extract just human alignment or simply base counting, we can specify what to count, coverage, common base, frequencies, and more
#this module was intended to use in junction with output files from teatime package, however it can be used as standalone
maf_matrix = extract_maf.maf_io(coordinate_table=alignment_coordinate, maf = MAF_dir, separated_maf=True, count_arg='common',custom_id=True)
#if there are teatime output files, they can be used to filter alignment to be counted
e_value_table_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/e_value/THE1C.txt'
internal_id_table_filepath = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/output/teatime/internal_id/THE1C.txt'
maf_matrix_filtered = extract_maf.maf_io(coordinate_table=alignment_coordinate, maf = MAF_dir, separated_maf=True, count_arg='common',custom_id=True, e_value_table=e_value_table_filepath, internal_id_table=internal_id_table_filepath)
#%% VCF file
# extract_vcf is a wrapper for cyvcf2
# as of now it supports gnomAD vcf file. As gnomAD VCF store subset of data with key such as AF (allele freqeucy), we can specify query key for a subset we want
from ma_mapper import extract_vcf
VCF_dir = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/allele_frequency_vcf/gnomad'
vcf_matrix = extract_vcf.vcf_io(coordinate_table=alignment_coordinate, vcf=VCF_dir, query_key='AF', vcf_format='gnomad')
#%%