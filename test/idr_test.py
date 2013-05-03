import unittest
from spp_idr import idr
import os
import shutil
import glob
from cluster_helper.cluster import cluster_view

dirname = os.path.dirname(__file__)
data_dir = os.path.join(dirname, "data")
test_dir = os.path.join(dirname, "test_dir")
peak_dir = os.path.join(data_dir, "peaks")
spp_path = os.path.join(dirname, "../scripts/phantompeakqualtools/run_spp.R")
idr_path = os.path.join(dirname, "../scripts/idrCode/")
idr_runner_path = os.path.join(idr_path, "batch-consistency-analysis.r")
idr_plotter_path = os.path.join(idr_path, "batch-consistency-plot.r")
control_file = os.path.join(data_dir, "control.tagAlign")
exp_file = os.path.join(data_dir, "ctcf.tagAlign")
peak_files = [os.path.join(peak_dir, "ctcf_0_0_VS_control_0_pooled.regionPeak.gz"),
              os.path.join(peak_dir, "ctcf_0_1_VS_control_0_pooled.regionPeak.gz")]
pooled_file = os.path.join(peak_dir, "ctcf_0_0_VS_control_0_pooled.regionPeak.gz")
scheduler = "lsf"
queue = "hsph"
jobs = 5
NPEAKS = 300000


class TestBamfileOperations(unittest.TestCase):
    def setUp(self):
        idr.safe_makedir(test_dir)
        bam_file = os.path.join(data_dir, "control.bam")
        self.bam_file = _copy_file_to_testdir(bam_file)

    def test_is_bamfile(self):
        self.assertTrue(idr.is_bamfile(self.bam_file))

    def test_bam_to_tagalign(self):
        tagalign_file = idr.bam_to_tagalign(self.bam_file)
        self.assertTrue(idr.is_tagalign(tagalign_file))
        self.assertTrue(idr.file_exists(tagalign_file))

        #    def tearDown(self):
        #shutil.rmtree(test_dir)


class TestTagAlignOperations(unittest.TestCase):

    def setUp(self):
        split_files = [control_file, exp_file]
        self.split_files = map(_copy_file_to_testdir, split_files)
        self.tagalign_file = split_files[0]

    def test_is_tagalign(self):
        self.assertTrue(idr.is_tagalign(self.tagalign_file))

    def test_tagalign_split(self):
        tagalign_split_files = idr.tagalign_split(self.tagalign_file)
        self.assertTrue(all(map(idr.is_tagalign, tagalign_split_files)))
        self.assertTrue(all(map(idr.file_exists, tagalign_split_files)))

    def test_tagalign_pool(self):
        pooled_file = idr.tagalign_pool(self.split_files)
        self.assertTrue(idr.is_tagalign(pooled_file))
        self.assertTrue(idr.file_exists(pooled_file))

        #    def tearDown(self):
        #shutil.rmtree(test_dir)

class TestSPPPeakCaller(unittest.TestCase):

    def setUp(self):
        idr.safe_makedir(test_dir)
        self.control = _copy_file_to_testdir(control_file)
        self.experimental = _copy_file_to_testdir(exp_file)

    def test_spp_caller(self):
        peak_caller = idr.spp_peak_caller(spp_path, cores=1)
        peaks = peak_caller(self.control, self.experimental)
        self.assertTrue(idr.file_exists(peaks))
        self.assertTrue(idr.is_peak_file(peaks))

        #    def tearDown(self):
        #shutil.rmtree(test_dir)

class TestCombinationCallers(unittest.TestCase):
    def setUp(self):
        self.control = _copy_file_to_testdir(control_file)
        self.experimental = _copy_file_to_testdir(exp_file)

    def test_pooled_caller(self):
        peak_caller = idr.spp_peak_caller(spp_path)
        control_replicates = idr.tagalign_split(self.control, nfiles=2)
        experimental_replicates = idr.tagalign_split(self.experimental,
                                                     nfiles=2)
        peaks = idr.call_peaks_on_pooled_replicates(control_replicates,
                                                    experimental_replicates,
                                                    peak_caller)
        self.assertTrue(idr.file_exists(peaks))
        self.assertTrue(idr.is_peak_file(peaks))

    def test_individual_caller(self):
        peak_caller = idr.spp_peak_caller(spp_path)
        control_replicates = idr.tagalign_split(self.control, nfiles=2)
        experimental_replicates = idr.tagalign_split(self.experimental,
                                                     nfiles=2)
        with cluster_view(scheduler, queue, jobs) as view:
            peaks = idr.call_peaks_on_individual_replicates(control_replicates,
                                                            experimental_replicates,
                                                            peak_caller, view.map)
        self.assertTrue(all(map(idr.file_exists, peaks)))
        self.assertTrue(all(map(idr.is_peak_file, peaks)))

    def test_pseudo_replicate_caller(self):
        peak_caller = idr.spp_peak_caller(spp_path)
        control_replicates = idr.tagalign_split(self.control, nfiles=2)
        experimental_replicates = idr.tagalign_split(self.experimental,
                                                     nfiles=2)
        #with cluster_view(scheduler, queue, jobs) as view:
        peaks = idr.call_peaks_on_pseudo_replicates(control_replicates,
                                                        experimental_replicates,
                                                        peak_caller, map)
        for replicate in peaks:
            self.assertTrue(all(map(idr.file_exists, replicate)))
            self.assertTrue(all(map(idr.is_peak_file, replicate)))

    def test_pooled_pseudo_replicate_caller(self):
        peak_caller = idr.spp_peak_caller(spp_path)
        control_replicates = idr.tagalign_split(self.control, nfiles=2)
        experimental_replicates = idr.tagalign_split(self.experimental,
                                                     nfiles=2)
        with cluster_view(scheduler, queue, jobs) as view:
            peaks = idr.call_peaks_on_pooled_pseudoreplicates(control_replicates,
                                                            experimental_replicates,
                                                            peak_caller, view.map)
        self.assertTrue(all(map(idr.file_exists, peaks)))
        self.assertTrue(all(map(idr.is_peak_file, peaks)))

        #    def tearDown(self):
        #shutil.rmtree(test_dir)

def setUp(self):
    idr.safe_makedir(test_dir)

def test_run_idr():
    individual_peaks = glob.glob(os.path.join(peak_dir, 'ctcf_[0-9]_VS_*.gz'))
    individual_peaks = map(_copy_file_to_testdir, individual_peaks)
    pseudo_peaks = glob.glob(os.path.join(peak_dir, 'ctcf_[0-9]_[0-9]_VS_*.gz'))
    pseudo_peaks = map(_copy_file_to_testdir, pseudo_peaks)
    pseudo_peaks.sort()
    pseudo_peaks = idr.bunch(pseudo_peaks)
    pooled_pseudo = glob.glob(os.path.join(peak_dir, 'ctcf_0_pooled_[0-9]_VS_*.gz'))
    pooled_pseudo = map(_copy_file_to_testdir, pooled_pseudo)
    idr_runner = idr.idr_runner(idr_runner_path)
    indiv, pseudo, pool_pseudo = idr.run_idr(individual_peaks, pseudo_peaks,
                                             pooled_pseudo, idr_runner)

    print indiv
    print pseudo
    print pool_pseudo
    assert all(map(idr.file_exists, indiv))
    assert all(map(idr.file_exists, list(idr.flatten(pooled_pseudo))))
    assert all(map(idr.file_exists, pool_pseudo))

def test_idr_runner():
    tmp_peak_files = map(_copy_file_to_testdir, peak_files)
    idr_runner = idr.idr_runner(idr_runner_path)
    idr_out = idr_runner(tmp_peak_files)
    assert idr.file_exists(idr_out)

def test_idr_plotter():
    prefixes = _prepare_idr_files(os.path.join(peak_dir, "idr", "reps"))
    idr_plotter = idr.idr_plotter(idr_plotter_path)
    plot_out = idr_plotter(prefixes)
    assert idr.file_exists(plot_out)

def _prepare_idr_files(dir):
    idr_files = glob.glob(os.path.join(dir, "*"))
    sub_dir = os.path.basename(dir)
    idr_files = map(_copy_file_to_testdir, idr_files, [sub_dir] * len(idr_files))
    idr_overlap = [x for x in idr_files if "overlapped" in x]
    prefixes = [x.split("-overlapped")[0] for x in idr_overlap]
    return prefixes

def _prepare_idr_overlapped(dir):
    prefixes = _prepare_idr_files(dir)
    return [x + "-overlapped-peaks.txt" for x in prefixes]

def _get_idr_set():
    i_files = _prepare_idr_files(os.path.join(peak_dir, "idr", "reps"))
    p_files = _prepare_idr_files(os.path.join(peak_dir, "idr", "pseudo"))
    pp_files = _prepare_idr_files(os.path.join(peak_dir, "idr", "pooled_pseudo"))
    return i_files, p_files, pp_files

def _get_idr_overlapped_set():
    i_files = _prepare_idr_overlapped(os.path.join(peak_dir, "idr", "reps"))
    p_files = _prepare_idr_overlapped(os.path.join(peak_dir, "idr", "pseudo"))
    pp_files = _prepare_idr_overlapped(os.path.join(peak_dir, "idr", "pooled_pseudo"))
    return i_files, p_files, pp_files

def test_post_idr_plot():
    i_files, p_files, pp_files = _get_idr_set()
    idr_plotter = idr.idr_plotter(idr_plotter_path)
    out_files = idr.plot_idr_output(i_files, p_files, pp_files, idr_plotter)
    assert all(map(idr.file_exists, list(idr.flatten(out_files))))

def test_get_peak_count_thresholds():
    idr_set = _get_idr_overlapped_set()
    i_counts, p_counts, pp_counts = idr.get_peak_count_thresholds(idr_set, NPEAKS)
    print i_counts, p_counts, pp_counts

def test_filter_peaks():
    idr_set = _get_idr_overlapped_set()
    tmp_pooled_file = _copy_file_to_testdir(pooled_file)
    filtered_peaks = idr.filter_peaks(tmp_pooled_file, idr_set, NPEAKS)
    print filtered_peaks
    assert all(map(idr.file_exists, filtered_peaks))

def test_run_analysis_no_replicates():
    tmp_control = [_copy_file_to_testdir(control_file)]
    tmp_exp = [_copy_file_to_testdir(exp_file)]
    with cluster_view(scheduler, queue, jobs) as view:
        plots, filtered_files = idr.run_analysis(tmp_control, tmp_exp,
                                                 spp_path, idr_runner_path,
                                                 idr_plotter_path, view.map)

    all(map(idr.file_exists, filtered_files))

def test_run_analysis():
    tmp_control = idr.tagalign_split(_copy_file_to_testdir(control_file))
    tmp_exp = idr.tagalign_split(_copy_file_to_testdir(exp_file))
    with cluster_view(scheduler, queue, jobs) as view:
        plots, filtered_files = idr.run_analysis(tmp_control, tmp_exp,
                                                 spp_path, idr_runner_path,
                                                 idr_plotter_path, view.map)

    all(map(idr.file_exists, filtered_files))


def _copy_file_to_testdir(fname, sub_dir=""):
    out_dir = os.path.join(test_dir, sub_dir)
    idr.safe_makedir(out_dir)
    out_file = os.path.join(out_dir, os.path.basename(fname))
    shutil.copyfile(fname, out_file)
    return out_file


def generate_data_from_bamfiles():
    test_files = glob.glob(os.path.join(data_dir, "*.bam"))
    map(idr.bam_to_tagalign, test_files)


if __name__ == "__main__":
    generate_data_from_bamfiles()
