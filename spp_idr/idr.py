import os
import pysam
import sys
import random
import subprocess
import time
import itertools
import gzip
import sh

def flag_problem_replicates():
    pass

def count_pseudoreplicate_peaks(idr_files, npeaks):
    threshold = select_self_consistency_threshold(npeaks)
    counts = []
    for idr_set in idr_files:
        counts.append(map(count_peaks_passing_threshold, idr_set,
                          [threshold] * len(idr_set)))
    return counts

def count_replicate_peaks(idr_individual, npeaks):
    threshold = select_self_consistency_threshold(npeaks)
    return map(count_peaks_passing_threshold, idr_individual,
               [threshold] * len(idr_individual))

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
    counter = 0
    def line_passes_threshold(line):
        x = line.split("\t")
        print x
        return float(x[11]) > threshold

    with open(peak_file) as in_handle:
        for line in in_handle:
            if line_passes_threshold(line):
                counter += 1
    return counter

def run_idr(individual_peaks, pseudo_peaks, pooled_pseudo_peaks, idr_runner):
    idr_individual = idr_on_peak_files(individual_peaks, idr_runner, "reps")
    idr_pseudo = map(idr_on_peak_files, pseudo_peaks, [idr_runner] * len(pseudo_peaks),
                     ["pseudo"] * len(pseudo_peaks))
    idr_pooled_pseudo = idr_on_peak_files(pooled_pseudo_peaks,
                                          idr_runner,
                                          "pooled_pseudo")
    return idr_individual, idr_pseudo, idr_pooled_pseudo

def idr_on_peak_files(peak_files, idr_runner, out_subdir="reps"):
    peak_files = map(gunzip, peak_files)
    pairs_to_run = itertools.combinations_with_replacement(peak_files, 2)
    out_files = map(idr_runner, pairs_to_run, [out_subdir] * len(pairs_to_run))
    return out_files

def call_peaks(controls, experimental, peak_caller):
    individual_peaks = call_peaks_on_individual_replicates(controls,
                                                           experimental,
                                                           peak_caller)
    pooled_peaks = call_peaks_on_pooled_replicates(controls,
                                                   experimental,
                                                   peak_caller)
    pseudo_peaks = call_peaks_on_pseudo_replicates(controls,
                                                   experimental,
                                                   peak_caller)
    pooled_pseudo_peaks = call_peaks_on_pooled_pseudoreplicates(controls,
                                                                experimental,
                                                                peak_caller)
    return individual_peaks, pooled_peaks, pseudo_peaks, pooled_pseudo_peaks

def call_peaks_on_pseudo_replicates(controls, experimental, peak_caller):
    pooled_controls = tagalign_pool(controls)
    pseudo_replicates = map(tagalign_split, experimental)
    peaks = map(peak_caller, [pooled_controls] * len(pseudo_replicates),
                pseudo_replicates)
    return peaks

def call_peaks_on_individual_replicates(controls, experimental, peak_caller):
    pooled_controls = tagalign_pool(controls)
    peaks = map(peak_caller, [pooled_controls] * len(experimental), experimental)
    return peaks

def call_peaks_on_pooled_replicates(controls, experimental, peak_caller):
    pooled_controls = tagalign_pool(controls)
    pooled_experimental = tagalign_pool(experimental)
    peaks = peak_caller(pooled_controls, pooled_experimental)
    return peaks

def call_peaks_on_pooled_pseudoreplicates(controls, experimental, peak_caller):
    pooled_controls = tagalign_pool(controls)
    pooled_experimental = tagalign_pool(experimental)
    pseudo_replicates = tagalign_split(pooled_experimental)
    peaks = map(peak_caller, [pooled_controls] * len(pseudo_replicates),
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


def _get_subdir(in_file, dirname):
    in_file_dir = os.path.dirname(in_file)
    return os.path.join(in_file_dir, dirname)

def bam_to_tagalign(in_file, out_file=None):
    base, _ = os.path.splitext(in_file)
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

def gunzip(fname):
    base, ext = os.path.splitext(fname)
    if ext != ".gz":
        return fname
    else:
        with open(fname, 'rb') as in_handle, open(base, 'wb') as out_handle:
            out_handle.writelines(in_handle)
        return base

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

    def __init__(self, spp_path, npeaks, cores=1):
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
               "-x=-500:85 -p={cores}")
        cmd_string = cmd.format(**locals())
        return cmd_string

    def _get_comparison_name(self, control_file, experimental_file):
        cbase, _ = os.path.splitext(control_file)
        ebase, _ = os.path.splitext(experimental_file)
        comparisons = map(os.path.basename, [ebase, cbase])
        return "_VS_".join(comparisons)

    def _get_stats_file(self, control_file, experimental_file):
        comparison_name = self._get_comparison_name(control_file,
                                                    experimental_file)
        stats_dir = _get_subdir(experimental_file, "stats")
        return os.path.join(stats_dir, comparison_name + ".tab")

    def _get_peaks_file(self, control, experimental):
        odir = _get_subdir(control, "peaks")
        out_file = (self._get_comparison_name(control, experimental) +
                    ".regionPeaks.gz")
        return os.path.join(odir, out_file)

    def __call__(self, control, experimental):
        cmd = self._spp_peak_call_command(control, experimental)
        peaks_file = self._get_peaks_file(control, experimental)
        stats_file = self._get_stats_file(control, experimental)
        self.map(self._call_peaks, [cmd])
        return peaks_file

class IDRRunner(object):
    def __init__(self, idr_path, results_subdir="reps"):
        self.idr_path = idr_path
        self.map = map
        self.results_subdir = results_subdir

    def set_map_function(self, mapper):
        self.map = mapper

    def _run_idr(self, pair):
        pair = map(gunzip, pair)
        out_file = self._get_out_file(pair)
        spp_path = self.spp_path
        p0 = pair[0]
        p1 = pair[1]
        cmd = ("Rscript {spp_path} {p0} {p1} -1 {out_file} 0 F signal.value")
        cmd_string = cmd.format(**locals())
        subprocess.check_call(cmd_string, shell=True)
        return out_file

    def _get_out_file(self, pair):
        parent_dir = os.path.dirname(os.path.dirname(pair[0]))
        out_dir = os.path.join(parent_dir, "idr", self.results_subdir)
        safe_makedir(out_dir)
        p1 = pair[0].split("_")[0]
        p2 = pair[1].split("_")[0]
        return os.path.join(out_dir, "_VS_".join([p1, p2]))

    def __call__(self, peak_files):
        idr_file = self.map(self._run_idr, [peak_files])
        return idr_file

def IDRPlotter(object):
    def __init__(self, idr_plotter_path):
        self.idr_plotter_path
        self.map = map

    def set_map_function(self, mapper):
        self.map = mapper

    def _run_idrplotter(self, idr_files, out_file):
        idr_plotter_path = self.idr_plotter_path
        idr_string = " ".join(idr_files)
        n = len(idr_files)
        cmd = ("Rscript {idr_plotter_path} {n} {out_file} {idr_string}")
        cmd_string = cmd.format(**locals())
        subprocess.check_call(cmd_string, shell=True)
        return out_file

    def __call__(self, idr_files, out_file):
        idr_plot = self.map(self._run_idrplotter, idr_files, [out_file])
        return idr_plot
