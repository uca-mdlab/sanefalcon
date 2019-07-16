import unittest
import os
import combine_profiles
import configparser
import shutil
import warnings


class TestCombineProfiles(unittest.TestCase):

    def setUp(self):
        self.upfiles = ['./tests/data/sample.bam.1.fwd', './tests/data/sample.bam.1.rev']
        self.downfiles = ['./tests/data/sample.bam.1.ifwd', './tests/data/sample.bam.1.irev']

    def test_get_profile_files_per_chrom(self):
        upfiles = list(combine_profiles.get_profile_files_per_chrom(self.upfiles, stream='up'))
        downfiles = list(combine_profiles.get_profile_files_per_chrom(self.downfiles, stream='down'))
        self.assertEqual(len(upfiles), len(downfiles))
        self.assertEqual(len(upfiles), 22)
        chrom1 = list(filter(lambda x: x[0] == 1, upfiles))
        self.assertEqual(len(chrom1), 1)
        self.assertEqual(len(chrom1[0]), 3)
        self.assertEqual(len(chrom1[0][1]), 1)

    def test_get_profile_file(self):
        res = combine_profiles.get_profile_file(self.upfiles[0])
        self.assertEqual(len(res), 147)
        self.assertEqual(res[0], 1981)

    def test_create_upstream(self):
        res = combine_profiles.create_upstream(self.upfiles)
        self.assertEqual(len(res['sample.bam']), 147)

    def test_create_dowstream(self):
        res = combine_profiles.create_downstream(self.downfiles)
        self.assertEqual(len(res['sample.bam']), 147)


if __name__ == '__main__':
    unittest.main()
