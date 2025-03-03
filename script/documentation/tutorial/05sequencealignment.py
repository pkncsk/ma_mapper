#%% LOAD PACKAGE
from ma_mapper import sequence_alignment
#%% seqeuence alignment wrapper
# as explained in the main page, ma_mapper is a wrapper package designed to use other bioinformatic package to do a specific task which is mapping/overlaying genome-wide data on multiple alignement of TE. It is mostly compatible with input from external package. However, this package also offer some module to steamline/simplify task such as sequence alignment
# This package wrap MAFFT fucntion so if MAFFT is installed it can be used from inside ma_mapper.
#%% extract coordinates from repeatmasker tabel
# the input of sequence alignment function is a file of TE sequences in FASTA format. One of the simple ways to extract dna seqeunces in the human genome is to use TE coordinate in BED format to extract sequences with SeqIO.
# the user can use extract_coord_from_repeatmasker_table to extract TE coordinates from a repeatmasker table
coord_table=sequence_alignment.extract_coord_from_repeatmasker_table(
    subfamily = 'THE1C',
    repeatmasker_table = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/hg38_repeatlib2014/hg38.fa.out.tsv'
    )
#%% extract TE seqeunces from human genome sequence and save them into fasta format file
te_seqeunces=sequence_alignment.sequence_io(
    coordinate_table=coord_table,
    source_fasta='/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/human_genome_fasta/hg38_fasta/hg38.fa',
    save_to_file='/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/THE1C.fa')
# %% align TE sequences using MAFFT warpper
# mafft arguments like nthread,  nthreadtb, nthreadtit can be used, additional parameter can be used by additional command using mafft_arg= )
sequence_alignment.mafft_align(
    input_filepath='/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/THE1C.fa',
    nthread=6,
    output_filepath='/rds/project/rds-XrHDlpCeVDg/users/pakkanan/data/resource/repeatmasker_table/THE1C.align'
)
# the output alignment file (fasta format) is ready to use in the downstream analyses
# %%
