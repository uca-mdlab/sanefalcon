import unittest
import os
import nucleosome.get_profile_parallel as gpp


class TestGetProfile(unittest.TestCase):

    def setUp(self):
        self.chromfilefwd = './tests/data/chromfile.fwd'
        self.chromfilerev = './tests/data/chromfile.rev'
        self.nuclfile = './tests/data/nuclfile'
        self.configurations = {
            'fwd0': {'outfile': '{}.fwd'.format(self.chromfilefwd), 'rev': 0},
            'fwd1': {'outfile': '{}.ifwd'.format(self.chromfilefwd), 'rev': 1},
            'rev0': {'outfile': '{}.irev'.format(self.chromfilerev), 'rev': 0},
            'rev1': {'outfile': '{}.rev'.format(self.chromfilerev), 'rev': 1},
        }

    def test_cast_line_to_numbers(self):
        a = ['1', '2', '3.2', '4.3', '5']
        res = gpp.cast_line_to_numbers(a)
        self.assertEqual(res, [1, 2, 3.2, 4.3, 5])
        a = ['1', '2', '3.2', '4.3', 'a']
        with self.assertRaises(ValueError):
            res = gpp.cast_line_to_numbers(a)

    def test_load_data_fwd(self):
        lines, peaks_fwd, reads_fwd = gpp.load_data(self.nuclfile, self.chromfilefwd)
        self.assertEqual(len(lines), 173875)
        self.assertEqual(len(peaks_fwd), 161844)
        self.assertEqual(len(reads_fwd), 71999)

    def test_load_data_rev(self):
        lines, peaks_rev, reads_rev = gpp.load_data(self.nuclfile, self.chromfilerev)
        self.assertEqual(len(lines), 173875)
        self.assertEqual(len(peaks_rev), 161844)
        self.assertEqual(len(reads_rev), 74985)

    def test_process_fwd0(self):
        outfile = self.configurations['fwd0']['outfile']
        lines, peaks, reads = gpp.load_data(self.nuclfile, self.chromfilefwd)
        sumPeaks = gpp.process_forward(peaks, reads, outfile)
        self.assertEqual(len(sumPeaks), 147)
        self.assertEqual(sumPeaks[0], 270.0)
        self.assertEqual(sumPeaks[-1], 297.0)

    def test_process_fwd1(self):
        outfile = self.configurations['fwd1']['outfile']
        lines, peaks, reads = gpp.load_data(self.nuclfile, self.chromfilefwd)
        sumPeaks = gpp.process_reverse(peaks, reads, outfile)
        self.assertEqual(len(sumPeaks), 147)
        self.assertEqual(sumPeaks[0], 270.0)
        self.assertEqual(sumPeaks[-1], 363.0)

    def test_process_rev0(self):
        outfile = self.configurations['rev0']['outfile']
        lines, peaks, reads = gpp.load_data(self.nuclfile, self.chromfilerev)
        sumPeaks = gpp.process_forward(peaks, reads, outfile)
        self.assertEqual(len(sumPeaks), 147)
        self.assertEqual(sumPeaks[0], 264.0)
        self.assertEqual(sumPeaks[-1], 318.0)

    def test_process_rev1(self):
        outfile = self.configurations['rev1']['outfile']
        lines, peaks, reads = gpp.load_data(self.nuclfile, self.chromfilerev)
        sumPeaks = gpp.process_reverse(peaks, reads, outfile)
        self.assertEqual(len(sumPeaks), 147)
        self.assertEqual(sumPeaks[0], 259.0)
        self.assertEqual(sumPeaks[-1], 314.0)

    def tearDown(self):
        for key, conf in self.configurations.items():
            try:
                os.remove(conf['outfile'])
            except FileNotFoundError:
                pass


if __name__ == '__main__':
    unittest.main()