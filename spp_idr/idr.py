import os
import pysam
import sys
import random
import subprocess
import time
import itertools
import gzip
import collections
import csv
import pandas as pd

"""
From: http://genome.ucsc.edu/FAQ/FAQformat.html#format12
This format is used to provide called peaks of signal enrichment based on pooled,
normalized (interpreted) data. It is a BED6+4 format.

    chrom - Name of the chromosome (or contig, scaffold, etc.).
    chromStart - The starting position of the feature in the chromosome or scaffold.
        The first base in a chromosome is numbered 0.
    chromEnd - The ending position of the feature in the chromosome or scaffold.
        The chromEnd base is not included in the display of the feature. For example,
        the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100,
        and span the bases numbered 0-99.
    name - Name given to a region (preferably unique). Use '.' if no name is assigned.
    score - Indicates how dark the peak will be displayed in the browser (0-1000).
        If all scores were '0' when the data were submitted to the DCC, the DCC assigned
        scores 1-1000 based on signal value. Ideally the average signalValue per base
        spread is between 100-1000.
    strand - +/- to denote strand or orientation (whenever applicable). Use '.' if
        no orientation is assigned.
    signalValue - Measurement of overall (usually, average) enrichment for the region.
    pValue - Measurement of statistical significance (-log10). Use -1 if no pValue
        is assigned.
    qValue - Measurement of statistical significance using false discovery rate (-log10).
        Use -1 if no qValue is assigned.
    peak - Point-source called for this peak; 0-based offset from chromStart.
        Use -1 if no point-source called.

Here is an example of narrowPeak format:

track type=narrowPeak visibility=3 db=hg19 name="nPk" description="ENCODE narrowPeak Example"
browser position chr1:9356000-9365000
chr1    9356548 9356648 .       0       .       182     5.0945  -1  50
chr1    9358722 9358822 .       0       .       91      4.6052  -1  40
chr1    9361082 9361182 .       0       .       182     9.2103  -1  75
"""

NARROWPEAK_HEADER = ["chrom", "chromStart", "chromEnd", "name", "score",
                     "strand", "signalValue", "pValue", "qValue", "peak"]

IDR_HEADER = ["chr1", "start1", "stop1", "sig.value1", "chr2", "start2", "stop2",
              "sig.value2", "idr.local", "IDR"]

REPLICATE_PASS = "PASS"
REPLICATE_WARNING = "WARNING"
REPLICATE_FAIL = "FAIL"

REPLICATE_WARNING_THRESHOLD = 2
REPLICATE_FAIL_THRESHOLD = 20
NPEAKS = 300000

def run_analysis(control_files, experimental_files, spp_path,
                 idr_runner_path, idr_plotter_path, mapper):

    # convert to tagalign
    print "Converting BAM files to tagAlign if necessary."
    control, experimental = map_2d(bam_to_tagalign, [control_files,
                                                     experimental_files], mapper)

    print "Calling peaks, this will take a while."
    # call peaks
    peak_caller = spp_peak_caller(spp_path, cores=1)
    i_peaks, p_peaks, pseudo_peaks, pp_peaks = call_peaks(control,
                                                          experimental,
                                                          peak_caller,
                                                          mapper)
    # run idr
    print "Performing IDR analysis."
    idr_run = idr_runner(idr_runner_path)
    i_idr, pseudo_idr, pp_idr = run_idr(i_peaks, pseudo_peaks, pp_peaks, idr_run)
    idr_plot = idr_plotter(idr_plotter_path)
    plots = plot_idr_output(i_idr, pseudo_idr, pp_idr, idr_plot)

    print "Filtering peaks using the cutoffs determined by IDR."
    # filter peaks
    filtered_files = filter_peaks(p_peaks, (i_idr, pseudo_idr, pp_idr), peak_caller.npeaks)

    print "Fin."
    return plots, filtered_files


def call_peaks(controls, experimental, peak_caller, mapper=map):
    individual_peaks = call_peaks_on_individual_replicates(controls,
                                                           experimental,
                                                           peak_caller, mapper)
    pooled_peaks = call_peaks_on_pooled_replicates(controls,
                                                   experimental,
                                                   peak_caller, mapper)
    pseudo_peaks = call_peaks_on_pseudo_replicates(controls,
                                                   experimental,
                                                   peak_caller, mapper)
    pooled_pseudo_peaks = call_peaks_on_pooled_pseudoreplicates(controls,
                                                                experimental,
                                                                peak_caller, mapper)
    return individual_peaks, pooled_peaks, pseudo_peaks, pooled_pseudo_peaks


def filter_peaks(peak_file, idr_set, npeaks):
    peak_file = gunzip(peak_file[0])
    sorted_peak_file = sort_peak_file(peak_file)
    conservative, optimum = get_filter_thresholds(idr_set, npeaks)
    conservative_peak_file = filter_peak_file(sorted_peak_file, conservative,
                                              "-conservative")
    optimum_peak_file = filter_peak_file(sorted_peak_file, optimum, "-optimum")
    return conservative_peak_file, optimum_peak_file

def sort_peak_file(peak_file):
    base, ext = os.path.splitext(peak_file)
    out_file = base + "-sorted" + ext
    if file_exists(out_file):
        return out_file
    df = pd.io.parsers.read_csv(peak_file, sep="\t", names=NARROWPEAK_HEADER)
    df = df.sort("signalValue", ascending=False)
    df.to_csv(out_file, sep="\t", index=False, header=False)
    return out_file

def filter_peak_file(sorted_peak_file, threshold, suffix):
    base, ext = os.path.splitext(sorted_peak_file)
    out_file = base + suffix + ext
    if file_exists(out_file):
        return out_file
    df = pd.io.parsers.read_csv(sorted_peak_file, sep="\t", names=NARROWPEAK_HEADER)
    filtered = df.head(threshold)
    filtered.to_csv(out_file, sep="\t", index=False)
    return out_file

def get_filter_thresholds(idr_set, npeaks):
    replicate_idr, pseudoreplicate_idr, pooled_pseudo_idr = idr_set
    original_replicate_threshold  = _original_replicate_threshold(replicate_idr,
                                                                  npeaks)[0]
    pooled_pseudoreplicate_threshold = _pooled_pseudoreplicate_threshold(pooled_pseudo_idr,
                                                                         npeaks)[0]
    conservative = original_replicate_threshold
    optimum = max([original_replicate_threshold, pooled_pseudoreplicate_threshold])
    return conservative, optimum

def get_peak_count_thresholds(idr_set, npeaks):
    replicate_idr, pseudoreplicate_idr, pooled_pseudo_idr = idr_set
    original_replicate_threshold  = _original_replicate_threshold(replicate_idr,
                                                                  npeaks)
    consistency_thresholds = count_pseudoreplicate_peaks(pseudoreplicate_idr,
                                                         npeaks)
    pooled_pseudoreplicate_threshold = _pooled_pseudoreplicate_threshold(pooled_pseudo_idr,
                                                                         npeaks)
    return (original_replicate_threshold, consistency_thresholds,
            pooled_pseudoreplicate_threshold)


def report_problem_replicates(replicates_with_flags):
   for item in replicates_with_flags:
       if item[1] == REPLICATE_WARNING:
           sys.stderr.write("%s has a self-consistency out of line with the other "
                            "replicates.")
       elif item[1] == REPLICATE_FAIL:
            sys.stderr.write("%s has a self-consistency way out of line with the other "
                             "replicates.")

def _original_replicate_threshold(replicate_peak_files, npeaks):
    if not replicate_peak_files:
        return [0]
    peaks = count_replicate_peaks(replicate_peak_files, npeaks)
    return [max(peaks)]

def _pooled_pseudoreplicate_threshold(pooled_pseudo_peak_file, npeaks):
    return count_pooled_replicate_peaks(pooled_pseudo_peak_file, npeaks)


def flag_problem_replicates(idr_files, npeaks):
    peaks = count_pseudoreplicate_peaks(idr_files, npeaks)
    mean_peaks = float(sum(peaks)) / len(peaks)
    distances = [abs(x - mean_peaks) for x in mean_peaks]
    flags = map(_replicate_flag, distances)
    return zip(idr_files, flags)

def _replicate_flag(distance_from_mean):
    if distance_from_mean > REPLICATE_FAIL_THRESHOLD:
        return REPLICATE_FAIL
    elif distance_from_mean > REPLICATE_WARNING_THRESHOLD:
        return REPLICATE_WARNING
    else:
        return REPLICATE_PASS


def count_pseudoreplicate_peaks(idr_files, npeaks):
    threshold = select_self_consistency_threshold(npeaks)
    return map(count_peaks_passing_threshold, idr_files,
               [threshold] * len(idr_files))

def count_replicate_peaks(idr_individual, npeaks):
    threshold = select_self_consistency_threshold(npeaks)
    return map(count_peaks_passing_threshold, idr_individual,
               [threshold] * len(idr_individual))


def count_pooled_replicate_peaks(idr_pooled, npeaks):
    threshold = select_pooled_consistency_threshold(npeaks)
    return map(count_peaks_passing_threshold, idr_pooled,
               [threshold] * len(idr_pooled))

def select_pooled_consistency_threshold(npeaks, large_genome=True):
    if npeaks >= 150000 and large_genome:
        return 0.005
    else:
        return 0.01

def select_self_consistency_threshold(npeaks, large_genome=True, shallow=False):
    if shallow:
        return 0.1
    elif npeaks >= 150000 and large_genome:
        return 0.02
    elif npeaks < 150000 and large_genome:
        return 0.05
    else:
        return 0.02

def count_peaks_passing_threshold(peak_file, threshold):
    peak_file = gunzip(peak_file)
    counter = 0

    with open(peak_file) as in_handle:
        peakreader = csv.DictReader(in_handle, delimiter=" ",
                                    fieldnames=IDR_HEADER)
        # skip the header
        peakreader.next()
        for line in peakreader:
            if float(line['IDR']) < threshold:
                counter += 1
    return counter

def count_peaks(peak_file):
    peak_file = gunzip(peak_file)
    counter = 0
    with open(peak_file) as in_handle:
        peakreader = csv.DictReader(in_handle, delimiter="\t",
                                    fieldnames=NARROWPEAK_HEADER)
        for line in peakreader:
                counter += 1
    return counter

def run_idr(individual_peaks, pseudo_peaks, pooled_pseudo_peaks, idr_runner):
    idr_individual = idr_on_peak_files(individual_peaks, idr_runner, "reps")
    idr_pseudo = map(idr_on_peak_files, pseudo_peaks, [idr_runner] * len(list(pseudo_peaks)),
                     ["pseudo"] * len(list(pseudo_peaks)))
    idr_pseudo = list(flatten(idr_pseudo))
    idr_pooled_pseudo = idr_on_peak_files(pooled_pseudo_peaks,
                                          idr_runner,
                                          "pooled_pseudo")
    return idr_individual, idr_pseudo, idr_pooled_pseudo

def idr_on_peak_files(peak_files, idr_runner, subdir):
    peak_files = map(gunzip, peak_files)
    pairs = list(itertools.combinations_with_replacement(peak_files, 2))
    pairs_to_run = [x for x in pairs if x[0] != x[1]]
    out_files = map(idr_runner, pairs_to_run, [subdir] * len(pairs_to_run))
    return out_files

def call_peaks_on_pseudo_replicates(controls, experimental, peak_caller, mapper):
    pooled_controls = tagalign_pool(controls)
    pseudo_replicates = mapper(tagalign_split, experimental)
    flat_replicates = list(flatten(pseudo_replicates))
    peaks = mapper(peak_caller, [pooled_controls] * len(flat_replicates),
                flat_replicates)
    return bunch(peaks)

def call_peaks_on_individual_replicates(controls, experimental, peak_caller, mapper):
    pooled_controls = tagalign_pool(controls)
    peaks = mapper(peak_caller, [pooled_controls] * len(experimental), experimental)
    return peaks

def call_peaks_on_pooled_replicates(controls, experimental, peak_caller, mapper):
    pooled_controls = tagalign_pool(controls)
    pooled_experimental = tagalign_pool(experimental)
    peaks = mapper(peak_caller, [pooled_controls], [pooled_experimental])
    return peaks

def call_peaks_on_pooled_pseudoreplicates(controls, experimental, peak_caller, mapper):
    pooled_controls = tagalign_pool(controls)
    pooled_experimental = tagalign_pool(experimental)
    pseudo_replicates = tagalign_split(pooled_experimental)
    peaks = mapper(peak_caller, [pooled_controls] * len(pseudo_replicates),
                pseudo_replicates)
    return peaks

#def spp_peak_caller(spp_path, npeak=300000):
#    def peak_caller(control_file, experimental_file):
#        cmd = _spp_peak_call_command(control_file, experimental_file, spp_path,
#                                     npeak)
#        out_file = _get_comparison_name(control_file,
#                                        experimental_file) + ".tagAlign"
#        subprocess.check_call(cmd, shell=True)
#        return out_file
#    return peak_caller

# use this for now until IPython that can pickle closures is released
def spp_peak_caller(spp_path, npeak=300000, cores=1):
    peak_caller = SPPPeakCaller(spp_path, npeak, cores)
    return peak_caller

def idr_runner(idr_runner_path):
    runner = IDRRunner(idr_runner_path)
    return runner

def idr_plotter(idr_plotter_path):
    runner = IDRPlotter(idr_plotter_path)
    return runner

def plot_idr_output(individual, pseudo, pooled_pseudo, idr_plotter):
    indiv_plot = idr_plotter(individual)
    pseudo_plots = idr_plotter(pseudo)
    pooled_pseudo = idr_plotter(pooled_pseudo)
    return indiv_plot, pseudo_plots, pooled_pseudo

def _get_subdir(in_file, dirname):
    in_file_dir = os.path.dirname(in_file)
    return os.path.join(in_file_dir, dirname)

def bam_to_tagalign(in_file, out_file=None):
    base, ext = os.path.splitext(in_file)
    if ext == ".tagAlign":
        return in_file
    if out_file is None:
        out_file = base + ".tagAlign"

    if file_exists(out_file):
        return out_file

    with _open_sam_or_bam(in_file) as in_handle, open(out_file, "w") as out_handle:
        for read in in_handle:
            if read.is_unmapped:
                continue
            out_handle.write(_uniqueread_to_tagalign(read, in_handle))

    return out_file

def tagalign_split(in_file, nfiles=2):
    "split a tagalign file into nfiles number of files, randomly writing reads"
    out_files = _get_tagalign_split_files(in_file, nfiles)
    out_handles = [open(x, "w") for x in out_files]
    with open(in_file) as in_handle:
        for read in in_handle:
            out_handles[random.randint(0, nfiles - 1)].write(read)
    [x.close() for x in out_handles]
    return out_files

def tagalign_pool(in_files, out_file=None):
    if out_file is None:
        base, ext = os.path.splitext(in_files[0])
        out_file = base + "_pooled" + ext
    with open(out_file, "w") as out_handle:
        for f in in_files:
            with open(f) as in_handle:
                for read in in_handle:
                    out_handle.write(read)
    return out_file

def _get_tagalign_split_files(in_file, nfiles=2):
    base, ext = os.path.splitext(in_file)
    fmt = "{0}_{1}{2}"
    out_files = [fmt.format(base, x, ext) for x in range(nfiles)]
    return out_files


def _uniqueread_to_tagalign(read, in_handle):
    chrom = in_handle.getrname(read.tid)
    fields = [chrom, read.pos, read.pos + read.qlen, read.seq,
              1000, _get_strand_of_read(read)]
    fields = map(str, fields)
    return "\t".join(fields) + "\n"

def _get_strand_of_read(read):
    if read.is_reverse:
        return "-"
    else:
        return "+"

def is_peak_file(in_file):
    base, ext = os.path.splitext(in_file)
    if ext == ".gz":
        _, ext = os.path.splitext(base)
    return ext == ".regionPeak"

def is_tagalign(in_file):
    _, ext = os.path.splitext(in_file)
    return ext == ".tagAlign"

def is_bamfile(in_file):
    _, ext = os.path.splitext(in_file)
    return ext == ".bam"

def is_samfile(in_file):
    _, ext = os.path.splitext(in_file)
    return ext == ".sam"

def _open_sam_or_bam(in_file):
    if is_samfile(in_file):
        flag = "r"
    elif is_bamfile(in_file):
        flag = "rb"
    else:
        raise ValueError("in_file is not a samfile or a bamfile")
        sys.exit(1)
    return pysam.Samfile(in_file, flag)

def safe_makedir(dname):
    """Make a directory if it doesn't exist, handling concurrent race conditions.
    """
    if not dname:
        return dname
    num_tries = 0
    max_tries = 5
    while not os.path.exists(dname):
        # we could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(dname)
        except OSError:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(2)
    return dname

def file_exists(fname):
    """Check if a file exists and is non-empty.
    """
    return os.path.exists(fname) and os.path.getsize(fname) > 0

def gunzip(fname):
    if not is_gzip(fname):
        return fname
    else:
        base, _ = os.path.splitext(fname)
        with gzip.open(fname, 'rb') as in_handle, open(base, 'w') as out_handle:
            out_handle.writelines(in_handle)
        return base


def is_gzip(fname):
    _, ext = os.path.splitext(fname)
    if ext == ".gz":
        return True
    else:
        return False

# this is here for now to get around IPython < 1.0 not being able to pickle closures.
class SPPPeakCaller(object):

    def __init__(self, spp_path, npeaks=NPEAKS, cores=1):
        self.spp_path = spp_path
        self.npeaks = npeaks
        self.map = map
        self.cores = cores

    def set_map_function(self, mapper):
        self.map = mapper

    def _call_peaks(self, cmd):
        return subprocess.check_call(cmd, shell=True)

    def _spp_peak_call_command(self, control, experimental):
        odir = _get_subdir(control, "peaks")
        safe_makedir(odir)
        stats = self._get_stats_file(control, experimental)
        safe_makedir(os.path.dirname(stats))
        spp_path = self.spp_path
        npeak = self.npeaks
        cores = self.cores
        cmd = ("Rscript {spp_path} -c={experimental} -i={control} "
               "-npeak={npeak} -odir={odir} -savr -savp -rf -out={stats} "
               "-x=-500:85")
        # -p={cores}")
        cmd_string = cmd.format(**locals())
        return cmd_string

    def _get_comparison_name(self, control_file, experimental_file):
        cbase, _ = os.path.splitext(control_file)
        ebase, _ = os.path.splitext(experimental_file)
        comparisons = map(os.path.basename, [ebase, cbase])
        return "_VS_".join(comparisons)

    def _get_stats_file(self, control_file, experimental_file):
        sys.stderr.write(control_file)
        sys.stderr.write(str(experimental_file))
        comparison_name = self._get_comparison_name(control_file,
                                                    experimental_file)
        stats_dir = _get_subdir(experimental_file, "stats")
        return os.path.join(stats_dir, comparison_name + ".tab")

    def _get_peaks_file(self, control, experimental):
        odir = _get_subdir(control, "peaks")
        out_file = (self._get_comparison_name(control, experimental) +
                    ".regionPeak.gz")
        return os.path.join(odir, out_file)

    def __call__(self, control, experimental):
        cmd = self._spp_peak_call_command(control, experimental)
        peaks_file = self._get_peaks_file(control, experimental)
        if file_exists(peaks_file):
            return peaks_file
        #stats_file = self._get_stats_file(control, experimental)
        self.map(self._call_peaks, [cmd])
        return peaks_file

class IDRRunner(object):
    def __init__(self, idr_path):
        self.idr_path = idr_path
        self.idr_dir = os.path.dirname(idr_path)
        self.idr_file = os.path.basename(idr_path)
        self.map = map

    def set_map_function(self, mapper):
        self.map = mapper

    def _run_idr(self, pair, out_dir):
        pair = map(gunzip, pair)
        out_prefix = self._get_out_prefix(pair, out_dir)
        idr_path = self.idr_file
        p0 = pair[0]
        p1 = pair[1]
        cur_dir = os.getcwd()
        out_file = self._get_out_file(pair, out_dir)
        if file_exists(out_file):
            return out_file
        os.chdir(self.idr_dir)
        cmd = ("Rscript {idr_path} {p0} {p1} -1 {out_prefix} 0 F signal.value")
        cmd_string = cmd.format(**locals())
        subprocess.check_call(cmd_string, shell=True)
        os.chdir(cur_dir)
        return out_file

    def _get_out_prefix(self, pair, out_dir):
        base_dir = os.path.dirname(pair[0])
        out_dir = os.path.join(base_dir, "idr", out_dir)
        safe_makedir(out_dir)
        p1 = os.path.basename(pair[0]).split("_VS")[0]
        p2 = os.path.basename(pair[1]).split("_VS")[0]
        return os.path.join(out_dir, "_VS_".join([p1, p2]))

    def _get_out_file(self, pair, out_dir):
        prefix = self._get_out_prefix(pair, out_dir)
        return prefix + "-overlapped-peaks.txt"

    def __call__(self, peak_files, out_dir="reps"):
        idr_file = self._run_idr(peak_files, out_dir)
        return idr_file

class IDRPlotter(object):
    def __init__(self, idr_plotter_path):
        self.idr_plotter_path = idr_plotter_path
        self.idr_dir = os.path.dirname(idr_plotter_path)
        self.idr_plot_file = os.path.basename(idr_plotter_path)
        self.map = map

    def set_map_function(self, mapper):
        self.map = mapper

    def _get_out_prefix(self, idr_files):
        print "in plotter" + str(idr_files)
        base_dir = os.path.dirname(idr_files[0])
        pieces = [os.path.basename(x).split("-overlapped")[0] for x in idr_files]
        return os.path.join(base_dir, "_VS_".join(pieces))

    def _run_idrplotter(self, idr_files, prefix):
        idr_plotter_path = self.idr_plot_file
        idr_files = [x.split("-overlapped")[0] for x in idr_files]
        idr_string = " ".join(idr_files)
        n = len(idr_files)
        cur_dir = os.getcwd()
        out_file = prefix + "-plot.ps"
        if file_exists(out_file):
            return out_file
        os.chdir(self.idr_dir)
        cmd = ("Rscript {idr_plotter_path} {n} {prefix} {idr_string}")
        cmd_string = cmd.format(**locals())
        subprocess.check_call(cmd_string, shell=True)
        os.chdir(cur_dir)
        return out_file

    def __call__(self, idr_files):
        if not idr_files:
            return []
        prefix = self._get_out_prefix(idr_files)
        idr_plot = self._run_idrplotter(idr_files, prefix)
        return idr_plot

def flatten(l):
    """
    flatten an irregular list of lists
    example: flatten([[[1, 2, 3], [4, 5]], 6]) -> [1, 2, 3, 4, 5, 6]
    lifted from: http://stackoverflow.com/questions/2158395/

    """
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el,
                                                                   basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el

def bunch(iterable, n=2, fillvalue=None):
    "bunch(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return map(list, itertools.izip_longest(fillvalue=fillvalue, *args))


def download_to_dir(url, dirname, extract=True, remove=True):
    cur_dir = os.getcwd()
    safe_makedir(dirname)
    os.chdir(dirname)
    cl = ["wget", url]
    subprocess.check_call(cl)
    if extract:
        cl = ["tar", "-xzvpf", os.path.basename(url)]
        subprocess.check_call(cl)
    if remove:
        os.remove(os.path.basename(url))
    os.chdir(cur_dir)
    return os.path.basename(url)

def map_2d(fn, l, mapper=map):
    """
    maps a function over a list of lists and groups the results
    at the end in the form of the original list of lists
    x = [[1, 2], [3, 4, 5]]
    list(map_2d(float, x)) -> [[1.0, 2.0], [3.0, 4.0, 5.0]]
    """
    def get_size(item):
        if isinstance(item, basestring):
            return 1
        else:
            return len(item)
    orig_sizes = map(get_size, l)
    results = iter(mapper(fn, (list(flatten(l)))))
    grouped = []
    for s in orig_sizes:
        grouped.append(take(s, results))
    return grouped

def take(n, iterable):
    """Return first n items of the iterable as a list

        >>> take(3, range(10))
        [0, 1, 2]
        >>> take(5, range(3))
        [0, 1, 2]

    Effectively a short replacement for ``next`` based iterator consumption
    when you want more than one item, but less than the whole iterator.

    """
    return list(itertools.islice(iterable, n))
