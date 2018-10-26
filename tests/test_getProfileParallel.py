import unittest
import os
import getProfileParallel


class TestGetProfile(unittest.TestCase):

    def setUp(self):
        self.chromfilefwd = './tests/chromfile.fwd'
        self.chromfilerev = './tests/chromfile.rev'
        self.nuclfile = './tests/nuclfile'
        self.configurations = {
            'fwd0': {'outfile': '{}.fwd'.format(self.chromfilefwd), 'rev': 0},
            'fwd1': {'outfile': '{}.ifwd'.format(self.chromfilefwd), 'rev': 1},
            'rev0': {'outfile': '{}.irev'.format(self.chromfilerev), 'rev': 0},
            'rev1': {'outfile': '{}.rev'.format(self.chromfilerev), 'rev': 1},
        }
        self.peaks_fwd = []
        self.reads_fwd = []
        self.peaks_rev = []
        self.reads_rev = []

    def test_cast_line_to_numbers(self):
        a = ['1', '2', '3.2', '4.3', '5']
        res = getProfileParallel.cast_line_to_numbers(a)
        self.assertEqual(res, [1, 2, 3.2, 4.3, 5])
        a = ['1', '2', '3.2', '4.3', 'a']
        with self.assertRaises(ValueError):
            res = getProfileParallel.cast_line_to_numbers(a)

    def test_load_data_fwd(self):
        lines, peaks_fwd, reads_fwd = getProfileParallel.load_data(self.nuclfile, self.chromfilefwd)
        self.assertEqual(len(lines), 173875)
        self.assertEqual(len(peaks_fwd), 161844)
        self.assertEqual(len(reads_fwd), 71999)

    def test_load_data_rev(self):
        lines, peaks_rev, reads_rev = getProfileParallel.load_data(self.nuclfile, self.chromfilerev)
        self.assertEqual(len(lines), 173875)
        self.assertEqual(len(peaks_rev), 161844)
        self.assertEqual(len(reads_rev), 74985)

    def test_process_forward_fwd0(self):
        outfile = self.configurations['fwd0']['outfile']
        lines, peaks_fwd, reads_fwd = getProfileParallel.load_data(self.nuclfile, self.chromfilerev)
        sumPeaksFwd = getProfileParallel.process_forward(peaks_fwd, reads_fwd, outfile)
        self.assertEqual(len(sumPeaksFwd), 147)
        self.assertEqual(sumPeaksFwd[0], '270.0')

    # def test_fwd0(self):
    #     outfile = self.configurations['fwd0']['outfile']
    #     rev = self.configurations['fwd0']['rev']
    #     result = os.system("python3 {} {} {} {} {}".format(self.scriptname, self.nuclfile, self.chromfilefwd, rev, outfile))
    #     self.assertEqual(result, 0)
    #     with open(outfile, 'r') as resfile:
    #         arr = resfile.readlines()[0].split(',')
    #
    #     self.assertEqual(arr[0].strip(), '270.0')
    #     self.assertEqual(arr[-1].strip(), '297.0')
    #
    # def test_fwd1(self):
    #     outfile = self.configurations['fwd1']['outfile']
    #     rev = self.configurations['fwd1']['rev']
    #     result = os.system("python3 {} {} {} {} {}".format(self.scriptname, self.nuclfile, self.chromfilefwd, rev, outfile))
    #     self.assertEqual(result, 0)
    #     with open(outfile, 'r') as resfile:
    #         arr = resfile.readlines()[0].split(',')
    #
    #     self.assertEqual(arr[0].strip(), '270.0')
    #     self.assertEqual(arr[-1].strip(), '363.0')
    #
    # def test_rev0(self):
    #     outfile = self.configurations['rev0']['outfile']
    #     rev = self.configurations['rev0']['rev']
    #     result = os.system("python3 {} {} {} {} {}".format(self.scriptname, self.nuclfile, self.chromfilefwd, rev, outfile))
    #     self.assertEqual(result, 0)
    #     with open(outfile, 'r') as resfile:
    #         arr = resfile.readlines()[0].split(',')
    #
    #     self.assertEqual(arr[0].strip(), '270.0')
    #     self.assertEqual(arr[-1].strip(), '297.0')
    #
    # def test_rev1(self):
    #     outfile = self.configurations['rev1']['outfile']
    #     rev = self.configurations['rev1']['rev']
    #     result = os.system("python3 {} {} {} {} {}".format(self.scriptname, self.nuclfile, self.chromfilefwd, rev, outfile))
    #     self.assertEqual(result, 0)
    #     with open(outfile, 'r') as resfile:
    #         arr = resfile.readlines()[0].split(',')
    #
    #     self.assertEqual(arr[0].strip(), '270.0')
    #     self.assertEqual(arr[-1].strip(), '363.0')

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()