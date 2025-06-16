#function take parameters
#BAM, bigwig, or BED (segway), or VCF-based frequency of snps / indels
#put extra sequences on the side? 500 bp of context wouldn't hurt to see alignment quality
#can we run MAFFT locally to have a fully automated script?

#TODO - add summary plot at the top
#TODO - option to show surrounding, unaligned sequences
#TODO - show where we deleted columns, or add option to output the unfiltered
#TODO - fix parameters, such as filtering bed files based on score, zoom, vmin vmax, max sequences, MAFFT parameters
#TODO - separate bound vs non bound, or display on the side, as well as families
#TODO - auto merge files with the wireframe
#TODO - add a parameter for bed files - segway
#TODO - bigbed support?
#TODO - clustering of profiles, heatmaps
#TODO - test bam files
#TODO - some script to find enriched families
#need to accept a dict to translate / group some target_repeats
#ideally a 2 part process - give a bed so it can extract sequences / produce an alignment
#2nd step is accept an aligned file + metadata file, + a directory of signal to display

#TODO: check if we want to presort everytime or ask for sorted files
#TODO: save that bed file along the sequence file, reopen both as a parameter to skip the alignment
#TODO: update the way we find unbound, just get the vector of which one is bound and update the names of those? left outer join?
#maybe outputting a df with - presence of binding site, max signal, total signal, etc
#check for presence of an aligned file - take it as a parameter, maybe even return it

#function: plot wireframe
#TODO:add 'both' as a color option and determine file name ourselves

#output and save align_data and metadata etc
#calculate shift by cross correlation ourselves
#TODO: sorting by values, or clustering by values, maybe even UMAP by values
#calculate vmax automatically

import sys
#cleanup later what we need / don't need from here
sys.path.append('/home/pc575/CAMproject/')
import myfunctions

from Bio import SeqIO
import pybedtools
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import numpy as np
import pandas as pd

import pickle
import re
import fastcluster
from scipy.cluster.hierarchy import leaves_list
import random
import glob
import os
import pysam
import math
import pyBigWig
import gc

def seq2int(seq):
    converted_seq = np.array([base2int(x) for x in seq], dtype=np.uint8)
    return converted_seq

def base2int(base):
    if base == '-':
        return 0
    elif base == 'A':
        return 1
    elif base == 'C':
        return 2
    elif base == 'T':
        return 3
    elif base == 'G':
        return 4
    elif base != '-':
        print(base)
        return 5

def normal_array(width=1, sigma=1, odd=0):
    ''' Returns an array of the normal distribution of the specified width '''
    sigma2 = float(sigma) ** 2

    def normal_func(x):
        return math.exp(-x * x / (2 * sigma2))

     # width is the half of the distribution
    values = list(map(normal_func, range(-width, width + odd)))
    values = np.array(values, np.float)

    return values

def merge_target_bed_files(target_bed_files, labels=None):
    #TODO: try to do all of this without the weird temp files / dependencies
    #TODO: memory-efficient for large number of files / large files
    def change_name_bed(feature, name):
        feature.name = name
        return feature

    merge_file = "/home/pc575/CAMproject/temp/temp.bed"
    f = open(merge_file, "w")
    f.close()
    f = open(merge_file, "a")

    for i, bed in enumerate(target_bed_files):
        bedtool = pybedtools.BedTool(bed)
        name = labels[i]
        bedtool = bedtool.each(change_name_bed, name).saveas()
        for line in bedtool:
            if str(line):
                f.write(str(line))
        del bedtool

    #TODO: newline inserted somewhere for no reason, for now i just delete it manually
    f.close()

    merge_bedtool = pybedtools.BedTool(merge_file).remove_invalid().sort().saveas(merge_file)
    return merge_bedtool

def filter_regions(merge_bedtool, repeats, max_sequences=None, overlap=True, fraction_overlap=0.5, adjust_str='', verbose=False):
    def adjust_name(x, adjust_str=''):
        x[3] = x[3] + adjust_str
        return x
    if overlap is True:
        intersect = repeats.intersect(merge_bedtool, F=fraction_overlap, u=True, sorted=True).saveas()
    else:
        intersect = repeats.intersect(merge_bedtool, F=fraction_overlap, v=True, sorted=True).saveas()
    if verbose:
        print("Length intersect ", len(intersect))

    intersect = intersect.each(adjust_name, adjust_str=adjust_str).saveas()

    intersect = intersect.filter(lambda x: int(x[2]) - int(x[1]) > 100).saveas()
    if verbose:
        print("Removed short repeats - ", len(intersect), " sequences remaining")

    if max_sequences > 0:
        factor = int(math.floor(len(intersect) / max_sequences))
        if factor > 1:
            print("Item too long:", len(intersect), "cutting down to " + str(max_sequences) + " seqs")
            intersect = pybedtools.BedTool(intersect[::factor]).saveas()

            print("Max sequences ", len(intersect))
    return intersect

def new_filter_msa(temp_array, col_threshold = 0.1, col_content_threshold = 0.1, row_threshold = 0.1):
    #remove columns that have less than col_threshold
    col_counts = np.zeros(temp_array.shape[1])
    col_counts_max = np.zeros(temp_array.shape[1])
    row_ratio = np.zeros(temp_array.shape[0])
    for i, row in temp_array.iterrows():
        nonzeros = np.nonzero(row.values)[0]
        col_counts[nonzeros] += 1
        col_counts_max[nonzeros[0]:nonzeros[-1]] += 1
    col_filter = (col_counts > col_counts_max * col_threshold) & (col_counts >= temp_array.shape[0] * col_content_threshold)

    #need to calculate % of gaps between non-zeros
    for i, row in temp_array.loc[:, col_filter].iterrows():
        nonzeros = np.nonzero(row.values)[0]
        if nonzeros.size and float(nonzeros[-1] - nonzeros[0]) > 0:
            row_ratio[i] = len(nonzeros) / (float(nonzeros[-1] - nonzeros[0]))

    row_filter = row_ratio >= row_threshold

    col_counts = np.zeros(temp_array.shape[1])
    col_counts_max = np.zeros(temp_array.shape[1])

    for i, row in temp_array.loc[row_filter, :].iterrows():
        nonzeros = np.nonzero(row.values)[0]

        col_counts[nonzeros] += 1
        col_counts_max[nonzeros[0]:nonzeros[-1]] += 1
    col_filter = (col_counts > col_counts_max * col_threshold) & (col_counts >= temp_array.shape[0] * col_content_threshold)

    #remove rows that have lost more than row_threshold columns
    #clear remaining columns that have too little sequences after filtering for outliers

    return (row_filter, col_filter)

def MAFFT_align(target_sequence_file, target_output):
    #call MAFFT from commandline
    #TODO: check if target_output exists, skip if it does
    import subprocess
    #cmd = 'mafft --ep 0.123 --thread 7 --kimura 1 --reorder {seqs} > {out}'.format(seqs=target_sequence_file, out="/mnt/d/temp/temp-mafft-aligned.fasta")
    #TODO: detect thread number - 1
    #TODO: test https://mafft.cbrc.jp/alignment/software/tips.html
    cmd = ['mafft', '--thread', '11', '--reorder', '--memsave', '--quiet', target_sequence_file]
    try:
        output = subprocess.check_output(cmd)
        with open(target_output, 'wb') as file:
            file.write(output)
    except subprocess.CalledProcessError as e:
        print(e.cmd)
        print(e.output)
        print(e.returncode)
        print("Something went wrong in MAFFT")

def get_sequences(bedtool, target_sequence_file, genome):
    def fields2name(f):
        #chr:start-end(strand)#name
        f[3] = f[0] + ":" + f[1] + "-" + f[2] + "(" +f[5] +")" + "#" + f[3]
        return f

    #TODO: filename should be generated or be a parameter
    #TODO: not even sure we need the seqs variable assigned

    bedtool.sequence(fi=genome, s=True, name=True).save_seqs(target_sequence_file)
    #check if existing alignment file was passed - maybe even check earlier?

def parse_alignment(fasta_file):
    #TODO: roll own fasta parser to get rid of seqio dependency
    regex = re.compile('(.+?)::(.+):(.+)-(.+)\((.)\)')
    records = (r for r in SeqIO.parse(fasta_file, "fasta"))
    seq_count = 0
    for i, value in enumerate(records):
        if seq_count == 0:
            seq_length = len(value)
        seq_count = seq_count + 1


    records = (r for r in SeqIO.parse(fasta_file, "fasta"))

    #can put as a numpy array instead, no need for pandas anywhere
    metadata_df = pd.DataFrame(columns=['family', 'chrom', 'start', 'end', 'strand'])
    metadata_df["start"] = metadata_df["start"].astype("int")
    metadata_df["end"] = metadata_df["end"].astype("int")

    align_array = np.zeros((seq_count, seq_length), dtype=np.uint8)
    for i, value in enumerate(records):
        #regex that to get everything at once
        split_values = regex.split(value.description)
        align_array[i, :] = np.array(str(value.seq).upper(), 'c').view(np.uint8)
        metadata_df.loc[i] = split_values[1:-1]

    #reformat DNA nucleotides to a simple number scheme
    #gap = 45 = 0
    #A = 65 = 1
    #C = 67 = 2
    #T = 84 = 3
    #G = 71 = 4
    #anything that is non-ATCG gets assigned as 5
    align_array[align_array == 45] = 0
    align_array[align_array == 65] = 1
    align_array[align_array == 67] = 2
    align_array[align_array == 84] = 3
    align_array[align_array == 71] = 4
    align_array[align_array >= 5] = 5
    align_data_df = pd.DataFrame(align_array, dtype=np.uint8)

    return align_data_df, metadata_df

def filter_alignment(align_data, metadata, col_threshold=0.05, col_content_threshold=0.05, row_threshold=0.50, max_sequences=15000):
    print(align_data.shape[0], "sequences of length", align_data.shape[1], "before filtering")
    filters = new_filter_msa(align_data, col_threshold=col_threshold, col_content_threshold=col_content_threshold, row_threshold=row_threshold)
    #TODO: add a print of sequence length after filtering
    return filters

def sort_alignment_per_family(align_data, metadata):
    #lots of stuff to cleanup here as well
    #would be a lot simpler to keep the ordering from mafft, allow for reordering according to groups
    concat_df = pd.DataFrame()
    index = []
    cat_list = []
    for name, group in metadata.groupby('family'):
        print("Clustering", name, align_data.loc[group.index, :].shape)
        index += myfunctions.cluster_data(align_data.loc[group.index, :], row_metric='hamming', col_metric=None, LO_rows='MOLO-avg', row_method='average').index.tolist()
        cat_list.append(name)
        gc.collect()

    from pandas.api.types import CategoricalDtype

    align_data = align_data.reindex(index)
    metadata = metadata.reindex(index)

    metadata['sort'] = pd.Series(range(metadata.shape[0]), index=metadata.index)
    metadata = metadata[['chrom', 'start', 'end', 'family', 'sort', 'strand']]
    metadata["start"] = metadata["start"].astype("int")
    metadata["end"] = metadata["end"].astype("int")
    #metadata['family'] = metadata['family'].astype(CategoricalDtype(categories=cat_list, ordered=True))
    #align_data = align_data.reindex(metadata.sort_values(by=['family', 'sort order']).index.values)

    return align_data, metadata

def process_bed_signal_simple(metadata, signal):
    signal_bed = pybedtools.BedTool(signal)
    target_bed = pybedtools.BedTool.from_dataframe(metadata)
    metadata_result = target_bed.intersect(signal_bed, loj=True, sorted=True).to_dataframe().sort_values('score')

    metadata_result[["cat", "count"]] = metadata_result.apply(count_cat, axis=1, result_type='expand')

    return metadata_result.sort_values('count', ascending=False).drop_duplicates('score').sort_values('score')['cat'].values

def count_cat(p):
    return (p.blockSizes, min(p.itemRgb, p.end) - max(p.thickEnd, p.start))

def process_bed_signal(metadata, signal):
    signal_bed = pybedtools.BedTool(signal).cut([0,1,2,3,4,5])
    #leftjoin metadata with bed, for each row where we have an intersect, build an array of BED signal starting at 0 from metadata coordinates
    #doing it with unsorted for now, would be better if we could assume its sorted, or we have the original order
    #better design possible
    target_bed = pybedtools.BedTool.from_dataframe(metadata)
    #score is the old sort order here
    target_bed = target_bed.sort()
    metadata_result = target_bed.intersect(signal_bed, loj=True, sorted=True).to_dataframe().sort_values('score')
    data = []
    p = metadata_result.iloc[0]
    current = p.score
    window = np.zeros(p.end - p.start)
    for i, p in metadata_result.iterrows():
        if p.score != current:
            data.append(window)
            window = np.zeros(p.end - p.start)
            current = p.score
        if p.thickStart != '.':
            #readjust if outside bounds
            if p.thickEnd < p.start:
                p.thickEnd = p.start
            if p.itemRgb > p.end:
                p.itemRgb = p.end

            if p.strand == '+':
                start = p.thickEnd - p.start
                end = p.itemRgb - p.start
                window[start:end] += np.full(p.itemRgb - p.thickEnd, p.blockSizes, dtype=float)
            else:
                start = p.end - p.itemRgb
                end =  p.end - p.thickEnd
                window[start:end] += np.full(p.itemRgb - p.thickEnd, p.blockSizes, dtype=float)[::-1]

    data.append(window)

    return data

def plot_embedding(embedding, original_df=None, clusterer=None, name=None, target_krab="ZNF267", color=None, size=3, alpha=0.5):
    from bokeh.plotting import figure, show, output_file, ColumnDataSource, output_notebook, save
    from bokeh.models import LogColorMapper, LogTicker, ColorBar, HoverTool, BoxZoomTool
    from bokeh.palettes import Viridis256, RdBu
    from bokeh.transform import linear_cmap
    from bokeh.io import export_png
    import itertools
    import colorcet as cc

    TOOLS="hover,pan,wheel_zoom,box_zoom,undo,redo,reset,save"

    embedding_df = pd.DataFrame(embedding)
    embedding_df = embedding_df.set_index(original_df.T.index.values)

    offset = 2
    cluster_colors = [next(itertools.islice(itertools.cycle(cc.glasbey), x+offset, x+offset+1)) if x >= 0
                      else (0.2, 0.2, 0.2)
                      for x in clusterer.labels_]

    if alpha == 'cluster':
        alpha = [20.0 if x == 124
                      else 0.0
                      for x in clusterer.labels_]

    project_colors = [next(itertools.islice(itertools.cycle(cc.glasbey), x+offset, x+offset+1)) if x >= 0
                      else (0.2, 0.2, 0.2)
                      for x in temp.groupby("proj_accession_BioProject").ngroup()]

    source = ColumnDataSource(data=dict(
        x=embedding[:,0],
        y=embedding[:,1],
        desc=original_df.T.index.values,
        title=temp['Title'],
        sample_title=temp['Sample_title'],
        bioproject=temp['proj_accession_BioProject'],
        color=temp[target_krab + "_log2"],
        rpkm=temp[target_krab],
        rpkm_scaled=temp[target_krab + "_scaled"] * 15,
        log2_fold_project=np.nan_to_num(np.log2((temp[target_krab]+1) / (temp[target_krab+"_project_min"]+1))),
        log2_fold_project_abs=np.absolute(np.nan_to_num(np.log2((temp[target_krab]+1) / (temp[target_krab+"_project_min"]+1)))),
        radius=(np.clip(temp[target_krab], a_min=0, a_max=100)) / 6,
        cluster_colors=cluster_colors,
        project_colors=project_colors,
        cluster=clusterer.labels_,
   #     alpha=alpha,
    ))

    #mapper = linear_cmap(field_name='color', palette=Viridis256, low=0, high=10)
    #mapper = linear_cmap(field_name='radius', palette=Viridis256, low=0, high=10)
    mapper = linear_cmap(field_name='log2_fold_project', palette=Viridis256, low=0, high=6)

    TOOLTIPS = [
        ("name", "@desc"),
        ("title", "@title"),
        ("Sample title", "@sample_title"),
        ("Bioproject", "@bioproject"),
        ("Expression", "@rpkm"),
        ("Scaled Expression", "@rpkm_scaled"),
        ("Fold", "@log2_fold_project"),
        ("Cluster", "@cluster"),
    ]

    p = figure(plot_width=1700, plot_height=800, tools=TOOLS, tooltips=TOOLTIPS, active_drag='pan', active_scroll='wheel_zoom',
               title="UMAP of " + target_krab + " expression in Skymap - Avg: " + str(round(temp[target_krab].mean(), 3)) + " - Max: " + str(round(temp[target_krab].max(), 3)) + " RPKM", output_backend="webgl")

    p.circle('x', 'y', source=source,
             alpha=alpha,
             line_width=0.1,
             line_color='#666666',
             size=size,
             #size='alpha',
             color=color,
            )

    #color_bar = ColorBar(color_mapper=mapper['transform'], width=25, location=(0,0))
    color_bar = ColorBar(color_mapper=mapper['transform'], width=25, location=(0,0))

    p.add_layout(color_bar, 'right')
    #export_png(p, filename="/mnt/d/Virtual machine/Shared/Skymap/PNG/" + target_krab + ".png")
    #output_file("/mnt/d/Virtual machine/Shared/Skymap/Interactive/" + target_krab + ".html")
    if name is not None:
        output_file("/home/pc575/CAMproject/temp/" + name + ".html")
        save(p)
    else:
        show(p)
    del p
    del bioProjectAnnotDf
    del technical_meta_data_df
    del srs_rearranged
    del temp
    gc.collect()

def process_bigwig_signal(metadata, signal, summary=False):
    bw = pyBigWig.open(signal)

    if summary is False:
        data = []
        for i, p in metadata.iterrows():
            window = np.array(bw.values(p.chrom, int(p.start), int(p.end)))
            if p.strand == '+':
                data.append(window)
            else:
                data.append(window[::-1])
    else:
        data = metadata.apply(lambda p: bw.stats(p.chrom, p.start, p.end, exact=True)[0], axis=1, raw=True).values

    bw.close()
    return data

def process_bam_signal(metadata, signal, mode='smooth-min'):
    #TODO: update variable names
    options_width = 100
    options_smooth = 100
    options_shift = 5
    bam = pysam.AlignmentFile(signal, "rb")

    normal = normal_array(width=options_width, sigma=options_smooth, odd=1)
    data = []

    for i, p in metadata.iterrows():
        width = p.end - p.start
        profile_normals_forward = np.zeros(width + (2*options_width), dtype=np.float16)
        profile_normals_reverse = np.zeros(width + (2*options_width), dtype=np.float16)
        profile_reads_forward = np.zeros(width + (2*options_width), dtype=np.uint8)
        profile_reads_reverse = np.zeros(width + (2*options_width), dtype=np.uint8)

        window = bam.fetch(p.chrom, p.start, p.end)

        for almnt in window:
            if almnt.is_reverse is False:
                read_position = almnt.reference_start - p.start + options_shift + options_width
            else:
                read_position = almnt.reference_end - p.start - options_shift - 1 + options_width

            start_in_window = read_position - options_width
            end_in_window = read_position + options_width + 1

            if (start_in_window < 0) or (end_in_window > width + (2*options_width)):
                continue

            if almnt.is_reverse is False:
                profile_normals_forward[start_in_window:end_in_window] += normal
                profile_reads_forward[read_position] += 1
            else:
                profile_normals_reverse[start_in_window:end_in_window] += normal
                profile_reads_reverse[read_position] += 1

        if mode == 'smooth-min':
            signal = np.minimum(profile_normals_forward, profile_normals_reverse)

        elif mode == 'smooth-max':
            signal = np.maximum(profile_normals_forward, profile_normals_reverse)
        else:
            signal = np.maximum(profile_reads_forward, profile_reads_reverse)

        if p.strand == '-':
            signal = signal[::-1]
        signal = signal[options_width:options_width+width]
        #signal = np.divide(signal, signal.max())
        data.append(signal)

    return data

def process_signal_files(alignment, metadata, col_filter, signal_files, cmap=None, color_list=None, vmin=None, vmax=None, display=True, save_dir=None, figsize=None, mode=None):
    #accept a cmap, else generate random ones

    if color_list is None:
        color_list = ["Blues", "Greens", "Purples", "Reds", "Oranges"]
        factor = int((len(signal_files) / len(color_list)) + 1)
        color_list = color_list * factor

    for i, signal_file in enumerate(signal_files):
        print("Processing", signal_file)
        #detect file type
        signal_ext = signal_file.lower()
        if signal_ext.endswith(('.bw', '.bigwig')):
            data = process_bigwig_signal(metadata, signal_file)
        elif signal_ext.endswith(('.bam')):
            data = process_bam_signal(metadata, signal_file, mode)
        elif signal_ext.endswith(('.bed', '.bed.gz')):
            data = process_bed_signal(metadata, signal_file)
        elif signal_ext.endswith(('.vcf')):
            data = process_vcf_signal(metadata, signal_file)
        else:
            print("Format unrecognized for ", signal_file)
            data = None

        if data:
            if cmap is None:
                cmap = plt.cm.get_cmap(color_list[i])
                cmap.set_under("#FFFFFF")
            if save_dir:
                save_name = save_dir + os.path.basename(signal_file) + '.png'
            else:
                save_name = None
            plot_alignment_signal(data, alignment, metadata, col_filter, cmap=cmap, vmin=vmin, vmax=vmax, display=display, save_name=save_name, figsize=figsize)

def plot_alignment_wireframe(data, color=True, display=False, save_name=None, zoom=None, figsize=None):
    if color is True:
        cmap = colors.ListedColormap(['#CCCCCC', 'green', 'yellow', 'red', 'blue', 'black', 'black'])
    else:
        cmap = colors.ListedColormap(['#8D8D8D', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#FFFFFF', '#000000', '#000000'])
    bounds = [0, 1, 2, 3, 4, 5, 6]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    plot_alignment(data, cmap, norm=norm, zoom=zoom, display=display, save_name=save_name, figsize=figsize)

def plot_alignment_signal(data, alignment, metadata, col_filter, cmap, vmin=None, vmax=None, display=True, save_name=None, figsize=None):
    alignment = alignment.reindex(metadata.index.values)
    signal_combined = np.zeros(alignment.shape, np.float32)
    i = 0
    for index, row in alignment.iterrows():
        signal_combined[i, row > 0] = data[i]
        i += 1

    signal_combined = signal_combined[:, col_filter]
    plot_alignment(signal_combined, cmap, zoom=None, vmin=vmin, vmax=vmax, display=display, save_name=save_name, figsize=figsize)

def plot_alignment(data, cmap, norm=None, zoom=None, vmin=None, vmax=None, display=False, save_name=None, figsize=None):
    #re-use this bottom part
    plt.rcParams['axes.linewidth'] = 3
    plt.rcParams['axes.edgecolor'] = '#222222'

    data = np.nan_to_num(data)
    if figsize is None:
        figsize = tuple(ti/100 for ti in data.T.shape)

    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor('white')
    ax.imshow(data, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax, interpolation='none', aspect='auto')
    ax.set_ylabel('Sequences')

    ax.grid(False)
    ax.autoscale(False)
    plt.tight_layout()
    if save_name:
        plt.savefig(save_name, dpi=300, interpolation='none')
    if display is True:
        plt.show()
    plt.close(fig)

def main_function(repeats=None, genome=None, target_regions=None, signal_files=None, target_repeats=None, max_sequences=0):
    #need to accept a dict to translate / group some target_repeats
    #ideally a 2 part process - give a bed so it can extract sequences / produce an alignment
    #2nd step is accept an aligned file + metadata file, + a directory of signal to display
    if None in [repeats, genome, target_regions, signal_files, target_repeats]:
        print("Missing parameters - exiting")
        return

    #TODO: check if .saveas() is necessary, check if we want to presort everytime or ask for sorted files
    filtered_repeats_bedtool = pybedtools.BedTool(repeats).filter(lambda x: x if x[3] in target_repeats else None).saveas()

    #process target regions
    #put filename in name column
    labels = [os.path.basename(file) for file in target_regions]
    merge_bedtool = merge_target_bed_files(target_regions, labels)

    repeats_intersect_target = filter_regions(merge_bedtool, filtered_repeats_bedtool, max_sequences=max_sequences, overlap=True, adjust_str='&bound')

    #doubly useless - just get the vector of which one is bound and updatet the names of those? left outer join?
    #CHECK: need to reload merge_bedtool?
    repeats_no_intersect_target = filter_regions(merge_bedtool, filtered_repeats_bedtool, max_sequences=max_sequences, overlap=False, adjust_str='&unbound')

    #CHECK: pretty sure thats useless and filtered_repeats_bedtool is the same, without updated change in name
    #TODO: save that bed file along the sequence file, reopen both as a parameter to skip the alignment
    total_intersect = repeats_intersect_target.cat(repeats_no_intersect_target, postmerge=False).sort().saveas()
    pybedtools.cleanup()

    target_sequence_file = "/home/pc575/CAMproject/temp/temp-mafft.fasta"
    aligned_sequence_file = "/home/pc575/CAMproject/temp/temp-mafft-aligned.fasta"
    get_sequences(total_intersect, target_sequence_file)

    MAFFT_align(target_sequence_file, aligned_sequence_file)

    #maybe outputting a df with - presence of binding site, max signal, total signal, etc
    #check for presence of an aligned file - take it as a parameter, maybe even return it

    align_data, metadata = parse_alignment(aligned_sequence_file)

    align_data_filtered, metadata_filtered = filter_alignment(align_data, metadata)
    align_data_sorted, metadata_sorted = sort_alignment_per_family(align_data_filtered, metadata_filtered)

    #function: plot wireframe
    #TODO:add 'both' as a color option and determine file name ourselves
    plot_alignment_wireframe(align_data_sorted, color=False, display=True, save_name="/home/pc575/CAMproject/temp/test_wireframe_bw.png")
    plot_alignment_wireframe(align_data_sorted, color='True', display=True, save_name="/home/pc575/CAMproject/temp/test_wireframe_color.png")

    #process data files, detecting their types as we go
    process_signal_files(total_intersect, signal_files)

    #output and save align_data and metadata etc
    #calculate shift by cross correlation ourselves