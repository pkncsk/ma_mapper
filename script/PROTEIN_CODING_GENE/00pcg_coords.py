#%%
import sys
import pandas as pd
import numpy as np
sys.path.append('/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/')
from ma_mapper import sequence_alignment
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
source_fasta = '/home/pc575/rds/rds-kzfps-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/hg38_fasta/hg38.fa'
global records
records = SeqIO.to_dict(SeqIO.parse(open(source_fasta), 'fasta'))
#%%
main_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
gencode_filepath ='/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/gencodeV44_tbl.txt'
cds_filepath ='/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/gencodeV44_cds_tbl.txt'
gencode_tbl=pd.read_csv(gencode_filepath, sep='\t', low_memory=False)
cds_tbl=pd.read_csv(cds_filepath, sep='\t', low_memory=False)
knownCds=cds_tbl['#name']
gencode_filtered = gencode_tbl[gencode_tbl.chrom.isin(main_chr)]
gencode_cds = gencode_filtered[gencode_filtered['#name'].isin(knownCds)]
gencode_protcoding=gencode_cds[gencode_cds.proteinID.notna()]
# %%
set_length = 500 #front 500 back 500
test_length = []
seq_records = []
for idx in gencode_protcoding.index:
    gene_by_idx=gencode_protcoding[gencode_protcoding.index==idx]
    metadata_id=gene_by_idx['#name'].values[0]
    print(metadata_id)
    start_list=gene_by_idx.exonStarts.str.split(',').values[0]
    end_list=gene_by_idx.exonEnds.str.split(',').values[0]
    cdsStart=gene_by_idx.cdsStart.values[0]
    cdsEnd=gene_by_idx.cdsEnd.values[0]
    cds_start_list = []
    cds_end_list = []
    if cdsEnd == cdsStart:
        cdsEnd = gene_by_idx.txEnd.values[0]
        cdsStart = gene_by_idx.txStart.values[0]

    for i in range(len(start_list)):
        if start_list[i]:
            #remove out of bound first
            if (int(start_list[i])<cdsStart) & (int(end_list[i])<cdsStart):
                continue
            elif (cdsEnd<int(start_list[i])) & (cdsEnd<int(end_list[i])):
                continue
            elif int(start_list[i]) <= cdsStart < int(end_list[i]):
                cds_start_list.append(cdsStart)
                cds_end_list.append(int(end_list[i]))
            elif (int(start_list[i])>cdsStart) & (int(end_list[i])<=cdsEnd):
                cds_start_list.append(int(start_list[i]))
                cds_end_list.append(int(end_list[i]))
            elif int(start_list[i]) < cdsEnd <= int(end_list[i]):
                cds_start_list.append(int(start_list[i]))
                cds_end_list.append(int(cdsEnd))
                break
    seq_strings=[]
    for idx in range(len(cds_start_list)):
        chrom = gene_by_idx.chrom.values[0]
        start = cds_start_list[idx]
        end = cds_end_list[idx]
        strand = gene_by_idx.strand.values[0]
        seq_string = sequence_alignment.extract_sequence(chrom, start, end, strand, records)
        seq_strings.append(seq_string)
        seqname = '::'.join([metadata_id,chrom,str(cdsStart),str(cdsEnd),strand])
    if strand == '-':
        seq_strings.reverse()
    seq_record = SeqRecord(Seq(''.join(seq_strings)),seqname , '', '')
    seq_records.append(seq_record)

    #post-hoc check
    prot_coding_seq = ''.join(seq_strings)
    pcs_length = len(prot_coding_seq)
    if pcs_length >=1000:
        accumulate_length = 0
        front = []
        for i in range(len(cds_start_list)):
            frag_length = cds_end_list[i] - cds_start_list[i]
            if accumulate_length+frag_length < 500:
                accumulate_length = accumulate_length+frag_length
                front.append([chrom, cds_start_list[i], cds_end_list[i],strand,metadata_id,])
                #print(i, accumulate_length, cds_start_list[i], cds_end_list[i])
            else:
                leftover = 500 - accumulate_length
                front.append([chrom, cds_start_list[i], cds_start_list[i]+leftover,strand,metadata_id,])
                test_length.extend(front)
                #print(i, leftover, cds_start_list[i], cds_start_list[0]+leftover)
                break
        accumulate_length = 0
        back = []
        for i in reversed(range(len(cds_start_list))):
            frag_length = cds_end_list[i] - cds_start_list[i]
            if accumulate_length+frag_length < 500:
                accumulate_length = accumulate_length+frag_length
                back.append([chrom, cds_start_list[i], cds_end_list[i],strand,metadata_id,])
                #print(i, accumulate_length, cds_start_list[i], cds_end_list[i])
            else:
                leftover = 500 - accumulate_length
                back.append([chrom, cds_end_list[i]-leftover, cds_end_list[i],strand,metadata_id,])
                back.reverse()
                test_length.extend(back)
                #print(i, leftover, cds_end_list[i]-leftover, cds_end_list[i])
                break

output_filepath = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/alignment/pcg.fasta'
with open(output_filepath, "w") as output_handle:
    SeqIO.write(seq_records, output_handle, "fasta")   
testgene_coord = pd.DataFrame(test_length, columns=['chrom','start','end','strand','id'])
output_filepath = '/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker_4_0_5_repeatlib20140131/old_result_redo/coord_internal_id/pcg_sliced.txt'
testgene_coord.to_csv(output_filepath, index=False,sep='\t')
#%%
set_length = 500 #front 500 back 500
test_length = []
seq_records = []
gene_by_idx=gencode_protcoding[gencode_protcoding['#name']=='ENST00000699311.1']
metadata_id=gene_by_idx['#name'].values[0]
print(metadata_id)
start_list=gene_by_idx.exonStarts.str.split(',').values[0]
end_list=gene_by_idx.exonEnds.str.split(',').values[0]
cdsStart=gene_by_idx.cdsStart.values[0]
cdsEnd=gene_by_idx.cdsEnd.values[0]
cds_start_list = []
cds_end_list = []
if cdsEnd == cdsStart:
    cdsEnd = gene_by_idx.txEnd.values[0]
    cdsStart = gene_by_idx.txStart.values[0]

for i in range(len(start_list)):
    if start_list[i]:
        #remove out of bound first
        if (int(start_list[i])<=cdsStart) & (int(end_list[i])<=cdsStart):
            continue
        elif (cdsEnd<int(start_list[i])) & (cdsEnd<int(end_list[i])):
            continue
        elif int(start_list[i]) <= cdsStart < int(end_list[i]):
            cds_start_list.append(cdsStart)
            cds_end_list.append(int(end_list[i]))
        elif (int(start_list[i])> cdsStart) & (int(end_list[i])<=cdsEnd):
            cds_start_list.append(int(start_list[i]))
            cds_end_list.append(int(end_list[i]))
        elif int(start_list[i]) < cdsEnd <= int(end_list[i]):
            cds_start_list.append(int(start_list[i]))
            cds_end_list.append(int(cdsEnd))
            break
seq_strings=[]
for idx in range(len(cds_start_list)):
    chrom = gene_by_idx.chrom.values[0]
    start = cds_start_list[idx]
    end = cds_end_list[idx]
    strand = gene_by_idx.strand.values[0]
    seq_string = sequence_alignment.extract_sequence(chrom, start, end, strand, records)
    seq_strings.append(seq_string)
    seqname = '::'.join([metadata_id,chrom,str(cdsStart),str(cdsEnd),strand])
if strand == '-':
    seq_strings.reverse()
seq_record = SeqRecord(Seq(''.join(seq_strings)),seqname , '', '')
seq_records.append(seq_record)

#post-hoc check
prot_coding_seq = ''.join(seq_strings)
pcs_length = len(prot_coding_seq)
if pcs_length >=1000:
    accumulate_length = 0
    front = []
    for i in range(len(cds_start_list)):
        frag_length = cds_end_list[i] - cds_start_list[i]
        if accumulate_length+frag_length < 500:
            accumulate_length = accumulate_length+frag_length
            front.append([chrom, cds_start_list[i], cds_end_list[i],strand,metadata_id,])
            print(i, accumulate_length, cds_start_list[i], cds_end_list[i])
        else:
            leftover = 500 - accumulate_length
            front.append([chrom, cds_start_list[i], cds_start_list[i]+leftover,strand,metadata_id,])
            test_length.extend(front)
            print(i, leftover, cds_start_list[i], cds_start_list[i]+leftover)
            break
    accumulate_length = 0
    back = []
    for i in reversed(range(len(cds_start_list))):
        frag_length = cds_end_list[i] - cds_start_list[i]
        if accumulate_length+frag_length < 500:
            accumulate_length = accumulate_length+frag_length
            back.append([chrom, cds_start_list[i], cds_end_list[i],strand,metadata_id,])
            print(i, accumulate_length, cds_start_list[i], cds_end_list[i])
        else:
            leftover = 500 - accumulate_length
            back.append([chrom, cds_end_list[i]-leftover, cds_end_list[i],strand,metadata_id,])
            back.reverse()
            test_length.extend(back)
            print(i, leftover, cds_end_list[i]-leftover, cds_end_list[i])
            break
# %%
gencode_protcoding[gencode_protcoding['#name']=='ENST00000699311.1']
# %%
seq_strings=[]
for idx in range(len(test_length)):
    chrom = test_length[idx][0]
    start = test_length[idx][1]
    end = test_length[idx][2]
    strand = test_length[idx][3]
    seq_string = sequence_alignment.extract_sequence(chrom, start, end, strand, records)
    seq_strings.append(seq_string)
    seqname = '::'.join([metadata_id,chrom,str(cdsStart),str(cdsEnd),strand])
if strand == '-':
    seq_strings.reverse()
seq_record = SeqRecord(Seq(''.join(seq_strings)),seqname , '', '')
seq_records.append(seq_record)
# %%
