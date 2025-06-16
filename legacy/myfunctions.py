# import sys
# sys.path.append('/media/sf_Shared/IPythonNotebooks/common/')
# import myfunctions
from __future__ import division
import os
import sys

import fnmatch
import random
import traceback
import time
import multiprocessing
import signal
import pybedtools
import pandas as pd
import numpy as np
import threading
import queue
from io import StringIO
from IPython.display import clear_output
import shlex
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import leaves_list
import seaborn as sns
sns.set_context("poster")
import fastcluster
import matplotlib
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
from scipy.cluster.hierarchy import linkage, to_tree

try:
    from os import scandir, walk
except ImportError:
    from scandir import scandir, walk

import scipy


def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        print('Finished in %2.2f seconds' % \
              (te - ts))
        return result

    return timed

chr_dict = {"1": "chr1",
            "2": "chr2",
            "3": "chr3",
            "4": "chr4",
            "5": "chr5",
            "6": "chr6",
            "7": "chr7",
            "8": "chr8",
            "9": "chr9",
            "10": "chr10",
            "11": "chr11",
            "12": "chr12",
            "13": "chr13",
            "14": "chr14",
            "15": "chr15",
            "16": "chr16",
            "17": "chr17",
            "18": "chr18",
            "19": "chr19",
            "20": "chr20",
            "21": "chr21",
            "22": "chr22",
            "X": "chrX",
            "Y": "chrY",
            "MT": "chrM"}

chr_dict_rev = {"chr1": "1",
                "chr2": "2",
                "chr3": "3",
                "chr4": "4",
                "chr5": "5",
                "chr6": "6",
                "chr7": "7",
                "chr8": "8",
                "chr9": "9",
                "chr10": "10",
                "chr11": "11",
                "chr12": "12",
                "chr13": "13",
                "chr14": "14",
                "chr15": "15",
                "chr16": "16",
                "chr17": "17",
                "chr18": "18",
                "chr19": "19",
                "chr20": "20",
                "chr21": "21",
                "chr22": "22",
                "chrX": "X",
                "chrY": "Y",
                "chrM": "MT"}


def hg19_to_GRCh37(feature):
    if str(feature.chrom) in chr_dict_rev:
        feature.chrom = chr_dict_rev[str(feature.chrom)]
    elif str(feature.chrom) not in chr_dict:
        feature = ''
    return feature


def GRCh37_to_hg19(feature):
    if str(feature.chrom) in chr_dict:
        feature.chrom = chr_dict[str(feature.chrom)]
    elif str(feature.chrom) not in chr_dict_rev:
        feature = ''
    return feature


def find_files(dir, ext):
    matches = []
    for root, dirnames, filenames in walk(dir):
        for filename in fnmatch.filter(filenames, ext):
            filename = os.path.join(root, filename)
            if file_len(filename) > 0:
                matches.append(os.path.join(root, filename))
    return matches


def worker_multi(q, results_queue, set_condition, current_args1, current_args2, bedtools1, bedtools2, labels1, labels2, status, f, F, e, genome, precise_fractions):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    while not q.empty():
        q_args = q.get()
        i = q_args[0]
        j = q_args[1]
        A_fraction = 0
        B_fraction = 0
        ratio = 0
        hits_len = 0
        query_len = 0
        ref_len = 0
        p_value = 1

        set_condition.acquire()
        try:
            while i in current_args1 or j in current_args2:
                set_condition.wait()
            current_args1[i] = True
            current_args2[j] = True
        finally:
            set_condition.release()
        try:
            a = bedtools1[i]
            b = bedtools2[j]

            c = a.fisher(b, g=genome, f=f, F=F, e=e)
            # This is now fixed in pybedtools 0.7.4 - bugged in previous versions
            hits_len = c.table["in -a"]["in -b"]
            query_len = c.table["in -a"]["not in -b"] + hits_len
            ref_len = c.table["not in -a"]["in -b"] + hits_len
            table = [[c.table["in -a"]["in -b"], c.table["in -a"]["not in -b"]], [c.table["not in -a"]["in -b"], c.table["not in -a"]["not in -b"]]]
            test = scipy.stats.fisher_exact(table, alternative="two-sided")
            ratio = test[0]
            two_tail_p_value = test[1]

            test = scipy.stats.fisher_exact(table, alternative="less")
            left_tail_p_value = test[1]

            test = scipy.stats.fisher_exact(table, alternative="greater")
            right_tail_p_value = test[1]

            # ratio = c.ratio
            # two_tail_p_value = c.two_tail
            # right_tail_p_value = c.right_tail
            # left_tail_p_value = c.left_tail
            if query_len != 0:
                if precise_fractions is True:
                    hits_len = len(a.intersect(b, sorted=True, u=True, f=f, F=F, e=e))
                    A_fraction = hits_len / query_len
                else:
                    A_fraction = hits_len / query_len
            else:
                A_fraction = 0
            if ref_len != 0:
                if precise_fractions is True:
                    B_fraction = len(b.intersect(a, sorted=True, u=True, f=f, F=F, e=e)) / query_len
                else:
                    B_fraction = hits_len / ref_len
            else:
                B_fraction = 0
            # del a
            # del b
            # del c
            if hits_len > 0:
                results_queue.put((int(i), int(j), float(A_fraction), float(B_fraction), float(ratio), float(two_tail_p_value), float(right_tail_p_value), float(left_tail_p_value), int(hits_len), int(query_len), int(ref_len)))

        except:
            print(traceback.format_exc())
            print("Unexpected error:", sys.exc_info())
            print("Encountered an error analyzing", labels1[i], "and", labels2[j])

        set_condition.acquire()
        try:
            del current_args1[i]
            del current_args2[j]
            set_condition.notify_all()
        finally:
            set_condition.release()
        status.put((i, j))
        pybedtools.cleanup()
        q.task_done()


def change_name_bed(feature, name):
    feature.name = name
    return feature


def cluster_bed_pybedtools_multi(savefile, bed_files1, bed_files2, labels1, labels2, reset=False, numthreads=1, f=1e-9, F=1e-9, e=False, prefilter=False, genome="/media/sf_Shared/g1k.genome", precise_fractions=False):
    if reset is False:
        # reload from file
        print("Reloading from files")
        sys.stdout.flush()

        data1 = pd.Panel({"A length": pd.read_csv(savefile + "_a_length.csv", index_col=0),
                          "B length": pd.read_csv(savefile + "_b_length.csv", index_col=0),
                          "Intersect": pd.read_csv(savefile + "_intersect.csv", index_col=0),
                          "A fraction": pd.read_csv(savefile + "_a_fraction.csv", index_col=0),
                          "B fraction": pd.read_csv(savefile + "_b_fraction.csv", index_col=0),
                          "Enrichment ratio": pd.read_csv(savefile + "_enrichment.csv", index_col=0),
                          "Fisher two-tail p-value": pd.read_csv(savefile + "_two-tail-p-value.csv", index_col=0),
                          "Fisher right-tail p-value": pd.read_csv(savefile + "_right-tail-p-value.csv", index_col=0),
                          "Fisher left-tail p-value": pd.read_csv(savefile + "_left-tail-p-value.csv", index_col=0),
                          })
    else:
        print("Processing with parameters: f =", f, ", F =", F, ", e =", e, ", g =", genome)
        print(len(bed_files1), "items compared against", len(bed_files2), "items")
        sys.stdout.flush()

        my_manager = multiprocessing.Manager()
        results_queue = my_manager.Queue()
        set_condition = my_manager.Condition()
        current_args1 = my_manager.dict()
        current_args2 = my_manager.dict()
        bedtools1 = my_manager.dict()
        bedtools2 = my_manager.dict()
        q = my_manager.Queue()
        status = my_manager.Queue()

        items_list = set()

        if prefilter is True:
            print("Prefiltering")
            sys.stdout.flush()
            temp_bedtools1 = list()
            temp_bedtools2 = list()

            print("Preloading bedtools A")
            sys.stdout.flush()
            for i, bed in enumerate(bed_files1):
                bedtool = pybedtools.BedTool(bed)
                bedtool = bedtool.each(change_name_bed, str(i)).remove_invalid().saveas()
                temp_bedtools1.append(pybedtools.BedTool(bedtool))

            print("Preloading bedtools B")
            sys.stdout.flush()
            for i, bed in enumerate(bed_files2):
                bedtool = pybedtools.BedTool(bed)
                bedtool = bedtool.each(change_name_bed, str(i)).remove_invalid().saveas()
                temp_bedtools2.append(pybedtools.BedTool(bedtool))

            print("Merging bedtools")
            sys.stdout.flush()
            merge1 = temp_bedtools1[0].saveas()
            merge1 = merge1.cat(*temp_bedtools1[1:], postmerge=False).sort().saveas()
            merge2 = temp_bedtools2[0].saveas()
            merge2 = merge2.cat(*temp_bedtools2[1:], postmerge=False).sort().saveas()

            print("Intersecting")
            sys.stdout.flush()
            results = merge1.intersect(merge2, wa=True, wb=True, sorted=True, f=f, F=F, e=e).saveas()
            print("Intersecting done")
            sys.stdout.flush()
            for row in results:
                items_list.add((int(row[3]), int(row[9])))

            for row in temp_bedtools1:
                del row

            for row in temp_bedtools2:
                del row

            del merge1
            del merge2
            del results
            pybedtools.cleanup()
            print("Prefiltering done")
            sys.stdout.flush()

        else:
            # fill the queue
            for i, bed1 in enumerate(bed_files1):
                for j, bed2 in enumerate(bed_files2):
                    items_list.add((i, j))

        # we reload everything even if prefiltered (pybedtools streaming workaround)
        print("Preloading bedtools A")
        sys.stdout.flush()
        for i, bed in enumerate(bed_files1):
            bedtools1[i] = pybedtools.BedTool(bed)
        print("Preloading bedtools B")
        sys.stdout.flush()
        for i, bed in enumerate(bed_files2):
            bedtools2[i] = pybedtools.BedTool(bed)

        items_list = list(items_list)
        # add entropy to diminish the likelihood of 2 processes using the same 2 files at once
        random.shuffle(items_list)
        print("Filling queue")
        sys.stdout.flush()
        for item in items_list:
            q.put(item)

        perc = pyprind.ProgPercent(q.qsize())

        data1 = pd.Panel({"A length": pd.DataFrame(index=labels1, columns=labels2, dtype=np.longfloat).fillna(0),
                          "B length": pd.DataFrame(index=labels1, columns=labels2, dtype=np.longfloat).fillna(0),
                          "Intersect": pd.DataFrame(index=labels1, columns=labels2, dtype=np.longfloat).fillna(0),
                          "A fraction": pd.DataFrame(index=labels1, columns=labels2, dtype=np.longfloat).fillna(0),
                          "B fraction": pd.DataFrame(index=labels1, columns=labels2, dtype=np.longfloat).fillna(0),
                          "Enrichment ratio": pd.DataFrame(index=labels1, columns=labels2, dtype=np.longfloat).fillna(0),
                          "Fisher two-tail p-value": pd.DataFrame(index=labels1, columns=labels2, dtype=np.longfloat).fillna(1.0),
                          "Fisher right-tail p-value": pd.DataFrame(index=labels1, columns=labels2, dtype=np.longfloat).fillna(1.0),
                          "Fisher left-tail p-value": pd.DataFrame(index=labels1, columns=labels2, dtype=np.longfloat).fillna(1.0),
                          })

        print("Processing", len(items_list), "items")
        sys.stdout.flush()
        workers = []
        for i in range(numthreads):
            tmp = multiprocessing.Process(target=worker_multi, args=(q, results_queue, set_condition, current_args1, current_args2, bedtools1, bedtools2, labels1, labels2, status, f, F, e, genome, precise_fractions))
            tmp.start()
            workers.append(tmp)

        try:
            while any(i.is_alive() for i in workers):
                time.sleep(1)
                while not status.empty():
                    i, j = status.get()
                    # clear_output(wait=True)
                    perc.update(item_id=labels1[i] + " vs " + labels2[j])
                sys.stdout.flush()

            for worker in workers:
                worker.join()

            while not results_queue.empty():
                item = results_queue.get()
                i = item[0]
                j = item[1]
                data1["A fraction"].ix[i, j] = item[2]
                data1["B fraction"].ix[i, j] = item[3]
                data1["Enrichment ratio"].ix[i, j] = item[4]
                data1["Fisher two-tail p-value"].ix[i, j] = item[5]
                data1["Fisher right-tail p-value"].ix[i, j] = item[6]
                data1["Fisher left-tail p-value"].ix[i, j] = item[7]
                data1["Intersect"].ix[i, j] = item[8]
                data1["A length"].ix[i, j] = item[9]
                data1["B length"].ix[i, j] = item[10]

            data1["A length"].to_csv(savefile + "_a_length.csv")
            data1["B length"].to_csv(savefile + "_b_length.csv")
            data1["Intersect"].to_csv(savefile + "_intersect.csv")
            data1["A fraction"].to_csv(savefile + "_a_fraction.csv")
            data1["B fraction"].to_csv(savefile + "_b_fraction.csv")
            data1["Enrichment ratio"].to_csv(savefile + "_enrichment.csv")
            data1["Fisher two-tail p-value"].to_csv(savefile + "_two-tail-p-value.csv")
            data1["Fisher right-tail p-value"].to_csv(savefile + "_right-tail-p-value.csv")
            data1["Fisher left-tail p-value"].to_csv(savefile + "_left-tail-p-value.csv")


        except KeyboardInterrupt:
            for worker in workers:
                worker.terminate()
                worker.join()
        finally:
            for key in iter(bedtools1.keys()):
                del bedtools1[key]
            for key in iter(bedtools2.keys()):
                del bedtools2[key]
            my_manager.shutdown()
            pybedtools.cleanup(remove_all=True)
    return data1


def worker_intervalstats(q, results_queue, set_condition, current_args1, current_args2, bed1, bed2, domain, intersect, labels1, labels2, status, p):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    while not q.empty():
        q_args = q.get()
        i = q_args[0]
        j = q_args[1]

        set_condition.acquire()
        try:
            while i in current_args1 or j in current_args2:
                set_condition.wait()
            current_args1[i] = True
            current_args2[j] = True
        finally:
            set_condition.release()

        try:
            command1 = "/media/sf_NGS/IntervalStats/Intervalstats -q " + bed1[i] + " -r " + bed2[j] + " -d " + domain + " -o " + '/dev/stdout'
            exit, out, err = get_exitcode_stdout_stderr(command1)
            output_strings = err.split("\n")
            query_len = 0

            if output_strings[1].split(" ")[-1]:
                query_len = int(output_strings[1].split(" ")[-1])

            results1 = pd.read_table(StringIO(out), sep="\t", header=None, names=["Query interval", "Closest reference interval", "Length of query", "Distance", "Numerator", "Denominator", "p-value"])
            results1 = results1[results1["p-value"] <= p]
            if intersect is True:
                results1 = results1[results1["Distance"] < intersect_dist]
            if query_len != 0:
                results_queue.put((int(i), int(j), 100 * len(results1) / query_len))
        except:
            print(traceback.format_exc())
            print("Unexpected error:", sys.exc_info())
            print("Encountered an error analyzing", labels1[i], "and", labels2[j])

        set_condition.acquire()
        try:
            del current_args1[i]
            del current_args2[j]
            set_condition.notify_all()
        finally:
            set_condition.release()
        status.put((i, j))
        q.task_done()

def cluster_intervalstats(savefile, bed_files1, bed_files2, labels1, labels2, domain, p=0.01, numthreads=2, intersect=False, intersect_dist=0, reset=False, save_dir=None):
    if reset is False:
        # reload from file
        print("Reloading from files")
        sys.stdout.flush()
        data1=pd.read_csv(savefile + ".csv", index_col=0)
    else:
        print("Processing with parameters: p =", p, ", intersect =", intersect, ", intersect_dist =", intersect_dist, ", domain =", domain)
        print(len(bed_files1), "items compared against", len(bed_files2), "items")
        sys.stdout.flush()

        my_manager = multiprocessing.Manager()
        results_queue = my_manager.Queue()
        set_condition = my_manager.Condition()
        current_args1 = my_manager.dict()
        current_args2 = my_manager.dict()
        q = my_manager.Queue()
        status = my_manager.Queue()

        items_list = set()

        # fill the queue
        for i, bed1 in enumerate(bed_files1):
            for j, bed2 in enumerate(bed_files2):
                items_list.add((i, j))

        items_list = list(items_list)
        # add entropy to diminish the likelihood of 2 processes using the same 2 files at once
        random.shuffle(items_list)
        print("Filling queue")
        sys.stdout.flush()
        for item in items_list:
            q.put(item)

        perc = pyprind.ProgPercent(q.qsize())

        data1 = pd.DataFrame(index=labels1, columns=labels2, dtype=np.longfloat).fillna(0)

        print("Processing", len(items_list), "items")
        sys.stdout.flush()
        workers = []
        for i in range(numthreads):
            tmp = multiprocessing.Process(target=worker_intervalstats, args=(q, results_queue, set_condition, current_args1, current_args2, bed_files1, bed_files2, domain, intersect, labels1, labels2, status, p))
            tmp.start()
            workers.append(tmp)

        try:
            while any(i.is_alive() for i in workers):
                time.sleep(1)
                while not status.empty():
                    i, j = status.get()
                    # clear_output(wait=True)
                    perc.update(item_id=labels1[i] + " vs " + labels2[j])
                sys.stdout.flush()

            for worker in workers:
                worker.join()

            while not results_queue.empty():
                item = results_queue.get()
                i = item[0]
                j = item[1]
                data1.ix[i, j] = item[2]

            data1.to_csv(savefile)

        except KeyboardInterrupt:
            for worker in workers:
                worker.terminate()
                worker.join()
        finally:
            my_manager.shutdown()
    return data1

def cluster_bed(savefile, bed_files1, bed_files2, labels1, labels2, domain, p=0.01, numthreads=2, intersect=False, intersect_dist=0, reset=False, save_dir=None):
    if reset is False:
        # reload from file
        data1 = pd.read_csv(savefile + ".csv", index_col=0)
    else:
        q = queue.queue()
        status = queue.queue()
        filename = "/dev/stdout"
        # fill the queue
        for i, bed1 in enumerate(bed_files1):
            for j, bed2 in enumerate(bed_files2):
                # do we need to calculate twice?
                command1 = "/media/sf_NGS/IntervalStats/Intervalstats -q " + bed1 + " -r " + bed2 + " -d " + domain + " -o " + filename
                q.put([(i, j), command1])
        clear_output(wait=True)
        print(q.qsize(), "elements in queue")
        sys.stdout.flush()
        perc = pyprind.ProgPercent(q.qsize())

        def worker():
            while True:
                if q.empty():
                    return
                args = q.get()
                i = args[0][0]
                j = args[0][1]
                exit, out, err = get_exitcode_stdout_stderr(args[1])
                output_strings = err.split("\n")
                query_len = 0

                if output_strings[1].split(" ")[-1]:
                    query_len = int(output_strings[1].split(" ")[-1])

                results1 = pd.read_table(StringIO(out), sep="\t", header=None, names=["Query interval", "Closest reference interval", "Length of query", "Distance", "Numerator", "Denominator", "p-value"])
                results1 = results1[results1["p-value"] <= p]
                if intersect is True:
                    results1 = results1[results1["Distance"] < intersect_dist]
                if query_len != 0:
                    data1.ix[i, j] = 100 * len(results1) / query_len
                status.put((i, j))
                q.task_done()

        threads = [threading.Thread(target=worker) for _i in range(numthreads)]
        data1 = pd.DataFrame(index=labels1, columns=labels2, dtype=np.longfloat).fillna(0)
        for thread in threads:
            thread.start()
        while any(i.is_alive() for i in threads):
            time.sleep(0.1)
            while not status.empty():
                i, j = status.get()
                # clear_output(wait=True)
                perc.update(item_id=labels1[i] + " vs " + labels2[j])
                sys.stdout.flush()
        for thread in threads:
            thread.join()
        data1.to_csv(savefile + ".csv")
    # clear_output(wait=True)
    sys.stdout.flush()
    sys.stdout.flush()
    return data1


def get_exitcode_stdout_stderr(cmd):
    """
    Execute the external command and get its exitcode, stdout and stderr.
    """
    args = shlex.split(cmd)

    proc = Popen(args, stdout=PIPE, stderr=PIPE)
    out, err = proc.communicate()
    exitcode = proc.returncode

    return exitcode, out, err


def plot_clusters(data1, metric="euclidean", method="single", cluster="None", vmax=None, vmin=None, colormap=plt.cm.jet):
    data2 = data1.copy(deep=True).transpose()

    #cluster dimension 1
    neworder1 = range(len(data1.index.values))
    neworder2 = range(len(data1.columns.values))

    if cluster == "1":
        linkmat = fastcluster.linkage(data1, method=method, metric=metric, preserve_input=False)
        neworder1 = list(leaves_list(linkmat))
        del pmat
        del linkmat
    elif cluster == "2":
        linkmat = fastcluster.linkage(data2, method=method, metric=metric, preserve_input=False)
        neworder2 = list(leaves_list(linkmat))
        del pmat
        del linkmat
    elif cluster == "Both":
        linkmat = fastcluster.linkage(data1, method=method, metric=metric, preserve_input=False)
        neworder1 = list(leaves_list(linkmat))
        del pmat
        del linkmat
        linkmat = fastcluster.linkage(data2, method=method, metric=metric, preserve_input=False)
        neworder2 = list(leaves_list(linkmat))
        del pmat
        del linkmat

    data1 = data1.reindex(index=[data1.index.values[label] for label in neworder1], copy=False)
    data1 = data1.reindex(columns=[data1.columns.values[label] for label in neworder2], copy=False)

    fig, ax = plt.subplots()

    if vmax is None:
        vmax = data1.values.max()

    if vmin is None:
        vmin = data1.values.min()

    heatmap = ax.pcolormesh(data1.values, cmap=colormap, alpha=0.8, vmax=vmax, vmin=vmin)

    # fig = plt.gcf()
    fig.set_size_inches(min(40, max(len(neworder2) * 0.5, 10)), min(40, max(len(neworder1) * 0.5, 10)))
    # ax.set_frame_on(False)

    ax.set_yticks(np.arange(data1.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(data1.shape[1]) + 0.5, minor=False)

    # ax.invert_yaxis()
    ax.xaxis.tick_top()

    ax.set_xticklabels(data1.columns.values, minor=False)
    ax.set_yticklabels(data1.index.values, minor=False)
    # plt.tick_params(axis='both', which='major', labelsize=8)
    # plt.tick_params(axis='both', which='minor', labelsize=8)
    plt.colorbar(heatmap, label="% of ERE family bound", orientation="horizontal", spacing="proportional")
    plt.xticks(rotation=90)

    # ax.grid(False)
    ax.axis('tight')

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    plt.tight_layout()
    # plt.savefig('/media/sf_Shared/test.pdf', format='pdf')

    plt.show()


def beds_vs_bed(bed_files, bed, domain, p=0.01, intersect=False, intersect_dist=0, reset=False, savefile=None):
    # accept a bed list and one bed
    # compare all beds agains the one bed
    # build a dataframe with 0/1
    # return
    if reset is False and os.path.isfile(savefile):
        # reload from file
        data1 = pd.read_csv(savefile, index_col=0)
    else:
        q = queue.queue()
        filename = "/dev/stdout"
        # fill the queue
        for i, bed1 in enumerate(bed_files):
            command1 = "/media/sf_NGS/IntervalStats/Intervalstats -q " + bed1 + " -r " + bed + " -d " + domain + " -o " + filename
            q.put([command1, bed1])

        def worker():
            while True:
                if q.empty():
                    return
                args = q.get()
                exit, out, err = get_exitcode_stdout_stderr(args[0])
                # output_strings = err.split("\n")
                results1 = pd.read_table(StringIO(out), sep="\t", header=None, names=["Query interval", "Closest reference interval", "Length of query", "Distance", "Numerator", "Denominator", "p-value"])
                results_filtered = results1.drop_duplicates(subset="Closest reference interval")
                results_filtered.set_index(keys="Closest reference interval", inplace=True)
                name = os.path.splitext(os.path.basename(args[1]))[0].split("_")[0]
                results_filtered[name] = results_filtered["p-value"] <= p
                if intersect is True:
                    results_filtered[name] = (results_filtered[name]) & (results_filtered["Distance"] < intersect_dist)

                with lock:
                    # change to row indexes, col index = values
                    data1.loc[results_filtered.index, name] = results_filtered[name].astype(int)
                q.task_done()

        threads = [threading.Thread(target=worker) for _i in range(4)]
        lock = threading.Lock()
        bed_df = pd.read_table(bed, sep="\t", header=None)
        bed_df = bed_df.drop_duplicates(subset=[0, 1, 2])
        data1 = pd.DataFrame(index=bed_df.apply(lambda x: '%s:%s:%s' % (x[0], x[1], x[2]), axis=1).values, columns=[os.path.splitext(os.path.basename(bed1))[0].split("_")[0] for bed1 in bed_files])
        data1 = data1.fillna(0)
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        data1.to_csv(savefile)
    return data1


def process_bed_directory(bed_files, threshold=5.0, merge=False):
    blacklist = pybedtools.BedTool("/media/sf_Shared/bed/consensusBlacklist.bed")
    for bed in bed_files:
        b = pybedtools.BedTool(bed)
        b = b.intersect(blacklist, v=True).saveas()
        b = b.sort().saveas()
        b = b.filter(lambda x: float(x.fields[6]) >= threshold).saveas()
        if b.count() > 0 and merge is True:
            b = b.merge(d=1000).saveas()
        b = b.saveas(('.').join(bed.split('.')[:-1]) + "_filtered.bed")


def file_len(fname):
    statinfo = os.stat(fname)
    return statinfo.st_size


def line_count(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

parula_data = [[0.2081, 0.1663, 0.5292],
               [0.2116, 0.1898, 0.5777],
               [0.2123, 0.2138, 0.6270],
               [0.2081, 0.2386, 0.6771],
               [0.1959, 0.2645, 0.7279],
               [0.1707, 0.2919, 0.7792],
               [0.1253, 0.3242, 0.8303],
               [0.0591, 0.3598, 0.8683],
               [0.0117, 0.3875, 0.8820],
               [0.0060, 0.4086, 0.8828],
               [0.0165, 0.4266, 0.8786],
               [0.0329, 0.4430, 0.8720],
               [0.0498, 0.4586, 0.8641],
               [0.0629, 0.4737, 0.8554],
               [0.0723, 0.4887, 0.8467],
               [0.0779, 0.5040, 0.8384],
               [0.0793, 0.5200, 0.8312],
               [0.0749, 0.5375, 0.8263],
               [0.0641, 0.5570, 0.8240],
               [0.0488, 0.5772, 0.8228],
               [0.0343, 0.5966, 0.8199],
               [0.0265, 0.6137, 0.8135],
               [0.0239, 0.6287, 0.8038],
               [0.0231, 0.6418, 0.7913],
               [0.0228, 0.6535, 0.7768],
               [0.0267, 0.6642, 0.7607],
               [0.0384, 0.6743, 0.7436],
               [0.0590, 0.6838, 0.7254],
               [0.0843, 0.6928, 0.7062],
               [0.1133, 0.7015, 0.6859],
               [0.1453, 0.7098, 0.6646],
               [0.1801, 0.7177, 0.6424],
               [0.2178, 0.7250, 0.6193],
               [0.2586, 0.7317, 0.5954],
               [0.3022, 0.7376, 0.5712],
               [0.3482, 0.7424, 0.5473],
               [0.3953, 0.7459, 0.5244],
               [0.4420, 0.7481, 0.5033],
               [0.4871, 0.7491, 0.4840],
               [0.5300, 0.7491, 0.4661],
               [0.5709, 0.7485, 0.4494],
               [0.6099, 0.7473, 0.4337],
               [0.6473, 0.7456, 0.4188],
               [0.6834, 0.7435, 0.4044],
               [0.7184, 0.7411, 0.3905],
               [0.7525, 0.7384, 0.3768],
               [0.7858, 0.7356, 0.3633],
               [0.8185, 0.7327, 0.3498],
               [0.8507, 0.7299, 0.3360],
               [0.8824, 0.7274, 0.3217],
               [0.9139, 0.7258, 0.3063],
               [0.9450, 0.7261, 0.2886],
               [0.9739, 0.7314, 0.2666],
               [0.9938, 0.7455, 0.2403],
               [0.9990, 0.7653, 0.2164],
               [0.9955, 0.7861, 0.1967],
               [0.9880, 0.8066, 0.1794],
               [0.9789, 0.8271, 0.1633],
               [0.9697, 0.8481, 0.1475],
               [0.9626, 0.8705, 0.1309],
               [0.9589, 0.8949, 0.1132],
               [0.9598, 0.9218, 0.0948],
               [0.9661, 0.9514, 0.0755],
               [0.9763, 0.9831, 0.0538]]

parula_map = matplotlib.colors.LinearSegmentedColormap.from_list('parula', parula_data)

matplotlib.cm.register_cmap('parula', parula_map)

def cluster_and_optimal_leaf_ordering(data, method='single', metric='euclidean', LO_method=None):
    #Normal clustering, no optimal leaf ordering
    try:
        if LO_method == 'OLO':
            links = linkage(data.values, method=method, metric=metric, optimal_ordering=True)
            tree = to_tree(links)
        else:
            links = linkage(data.values, method=method, metric=metric, optimal_ordering=False)
            tree = to_tree(links)
            #Leaf ordering as 2001 paper
            #MOLO, 2015
            if LO_method == 'MOLO-min':
                tree = MOLO_sort_smallest(tree)
            elif LO_method == 'MOLO-avg':
                tree = MOLO_sort_average(tree)
        return tree.pre_order(lambda x: x.id)
    except:
        print("Clustering went wrong")

from copy import deepcopy
def MOLO_sort_smallest(d):
  left = d.left
  right = d.right

  if left.is_leaf() and right.is_leaf():
    return(d)
  elif ~left.is_leaf() and right.is_leaf():
    d.left = MOLO_sort_smallest(left)
    d.dist = min(d.dist, d.left.dist)
  elif left.is_leaf() and ~right.is_leaf():
    d.left = MOLO_sort_smallest(right)
    d.right = left
    d.dist = min(d.dist, d.left.dist)
  else:
    lft = MOLO_sort_smallest(left)
    rght = MOLO_sort_smallest(right)
    left.dist = lft.dist
    right.dist = rght.dist
    if (left.dist <= right.dist):
      d.left = lft
      d.right = rght
    else:
      d.left = rght
      d.right = lft

    d.dist = min(d.dist, left.dist, right.dist)
  return(d)

def MOLO_sort_average(d):
  left = d.left
  right = d.right

  if left.is_leaf() and right.is_leaf():
    return(d)
  elif ~left.is_leaf() and right.is_leaf():
    d.left = MOLO_sort_average(left)
    d.dist = d.dist + d.left.dist
  elif left.is_leaf() and ~right.is_leaf():
    d.left = MOLO_sort_average(right)
    d.right = left
    d.dist = d.dist + d.left.dist
  else:
    lft = MOLO_sort_average(left)
    rght = MOLO_sort_average(right)
    left.dist = lft.dist
    right.dist = rght.dist

    left_avg = left.dist / (left.count - 1)
    right_avg = right.dist / (right.count - 1)

    if (left_avg <= right_avg):
      d.left = lft
      d.right = rght
    else:
      d.left = rght
      d.right = lft

    d.dist = d.dist + left.dist + right.dist
  return(d)

def cluster_data(data, row_metric='euclidean', col_metric='euclidean', row_method='average', col_method='average', LO_rows=None, LO_cols=None, cluster_on='original', return_data_type='original'):
    target_data = data.copy(deep=True)
    return_data = target_data
    norm_data = target_data.sub(target_data.mean(axis=1), axis=0).div(target_data.std(axis=1), axis=0)
    if cluster_on == "normalized":
        target_data = norm_data
    if return_data_type == 'normalized':
        return_data = norm_data

    if col_metric is not None:
        reordered_cols = cluster_and_optimal_leaf_ordering(target_data, method=col_method, metric=col_metric, LO_method=LO_cols)
        return_data = return_data.iloc[reordered_cols, :]

    if row_metric is not None:
        reordered_rows = cluster_and_optimal_leaf_ordering(target_data.transpose(), method=row_method, metric=row_metric, LO_method=LO_rows)
        return_data = return_data.iloc[:, reordered_rows]

    return return_data
