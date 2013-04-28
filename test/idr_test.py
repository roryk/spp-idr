import unittest
from spp_idr import idr
import os
import shutil
import glob

dirname = os.path.dirname(__file__)
data_dir = os.path.join(dirname, "data")
test_dir = os.path.join(dirname, "test_dir")
spp_path = os.path.join(dirname, "../scripts/phantompeakqualtools/run_spp.R")
control_file = os.path.join(data_dir, "control.tagAlign")
exp_file = os.path.join(data_dir, "ctcf.tagAlign")


class TestBamfileOperations(unittest.TestCase):
    def setUp(self):
        idr.safe_makedir(test_dir)
        bam_file = os.path.join(data_dir, "control.bam")
        self.bam_file = self._copy_file_to_testdir(bam_file)

    def test_is_bamfile(self):
        self.assertTrue(idr.is_bamfile(self.bam_file))

    def test_bam_to_tagalign(self):
        tagalign_file = idr.bam_to_tagalign(self.bam_file)
        self.assertTrue(idr.is_tagalign(tagalign_file))
        self.assertTrue(idr.file_exists(tagalign_file))

    def _copy_file_to_testdir(self, fname):
        out_file = os.path.join(test_dir, os.path.basename(fname))
        shutil.copyfile(fname, out_file)
        return out_file

#    def tearDown(self):
#        shutil.rmtree(test_dir)


class TestTagAlignOperations(unittest.TestCase):

    def setUp(self):
        split_files = [control_file, exp_file]
        self.split_files = map(self._copy_file_to_testdir, split_files)
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

    def _copy_file_to_testdir(self, fname):
        out_file = os.path.join(test_dir, os.path.basename(fname))
        shutil.copyfile(fname, out_file)
        return out_file

#    def tearDown(self):
#        shutil.rmtree(test_dir)

class TestSPPPeakCaller(unittest.TestCase):

    def setUp(self):
        idr.safe_makedir(test_dir)
        self.control = self._copy_file_to_testdir(control_file)
        self.experimental = self._copy_file_to_testdir(exp_file)

    def test_spp_caller(self):
        peak_caller = idr.spp_peak_caller(spp_path, cores=10)
        peaks = peak_caller(self.control, self.experimental)
        self.assertTrue(idr.file_exists(peaks))
        self.assertTrue(idr.is_peak_file(peaks))

    def _copy_file_to_testdir(self, fname):
        idr.safe_makedir(test_dir)
        out_file = os.path.join(test_dir, os.path.basename(fname))
        shutil.copyfile(fname, out_file)
        return out_file

#    def tearDown(self):
#        shutil.rmtree(test_dir)

class TestPipeline(unittest.TestCase):
    def setUp(self):
        self.control = self._copy_file_to_testdir(control_file)
        self.experimental = self._copy_file_to_testdir(exp_file)

    def test_pooled_caller(self):
        peak_caller = idr.spp_peak_caller(spp_path)
        control_replicates = idr.tagalign_split(self.control, nfiles=5)
        experimental_replicates = idr.tagalign_split(self.experimental,
                                                     nfiles=5)
        peaks = idr.call_peaks_on_pooled_replicates(control_replicates,
                                                    experimental_replicates,
                                                    peak_caller)
        print peaks
        self.assertTrue(idr.file_exists(peaks))
        self.assertTrue(idr.is_peak_file(peaks))

    def _copy_file_to_testdir(self, fname):
        idr.safe_makedir(test_dir)
        out_file = os.path.join(test_dir, os.path.basename(fname))
        shutil.copyfile(fname, out_file)
        return out_file

    def test_individual_caller(self):
        peak_caller = idr.spp_peak_caller(spp_path)
        control_replicates = idr.tagalign_split(self.control, nfiles=5)
        experimental_replicates = idr.tagalign_split(self.experimental,
                                                     nfiles=5)
        peaks = idr.call_peaks_on_individual_replicates(control_replicates,
                                                        experimental_replicates,
                                                        peak_caller)
        print peaks
        self.assertTrue(all(map(idr.file_exists, peaks)))
        self.assertTrue(all(map(idr.is_peak_file, peaks)))

#    def test_pseudo_replicate_caller(self):
#        peak_caller = idr.spp_peak_caller(spp_path)
#        control_replicates = idr.tagalign_split(self.control, nfiles=5)
#        experimental_replicates = idr.tagalign_split(self.experimental,
#                                                     nfiles=5)
#        peaks = idr.call_peaks_on_pseudo_replicates(control_replicates,
#                                                    experimental_replicates,
#                                                    peak_caller)
#        print peaks
#        for replicate in peaks:
#            self.assertTrue(all(map(idr.file_exists, replicate)))
#            self.assertTrue(all(map(idr.is_peak_file, replicate)))
#
#    def test_pooled_pseudo_replicate_caller(self):
#        peak_caller = idr.spp_peak_caller(spp_path)
#        control_replicates = idr.tagalign_split(self.control, nfiles=5)
#        experimental_replicates = idr.tagalign_split(self,experimental,
#                                                     nfiles=5)
#        peaks = idr.call_peaks_on_pooled_pseudoreplicates(control_replicates,
#                                                          experimental_replicates,
#                                                          peak_caller)
#        print peaks
#        self.assertTrue(all(map(idr.file_exists, peaks)))
#        self.assertTrue(all(map(idr.is_peak_file, peaks)))

def generate_test_data_from_bamfiles():
    test_files = glob.glob(os.path.join(data_dir, "*.bam"))
    map(idr.bam_to_tagalign, test_files)

if __name__ == "__main__":
    generate_test_data_from_bamfiles()
