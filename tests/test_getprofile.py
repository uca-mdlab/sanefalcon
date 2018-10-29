import unittest
import os


class TestGetProfile(unittest.TestCase):
    def setUp(self):
        self.scriptname = 'getProfile.py'
        self.chromfilefwd = './tests/data/chromfile.fwd'
        self.chromfilerev = './tests/data/chromfile.rev'
        self.nuclfile = './tests/data/nuclfile'
        self.configurations = {
            'fwd0': {'outfile': '{}.p.fwd'.format(self.chromfilefwd), 'rev': 0},
            'fwd1': {'outfile': '{}.p.ifwd'.format(self.chromfilefwd), 'rev': 1},
            'rev0': {'outfile': '{}.p.irev'.format(self.chromfilerev), 'rev': 0},
            'rev1': {'outfile': '{}.p.rev'.format(self.chromfilerev), 'rev': 1},
        }

    def test_fwd0(self):
        outfile = self.configurations['fwd0']['outfile']
        rev = self.configurations['fwd0']['rev']
        result = os.system("python3 {} {} {} {} {}".format(self.scriptname, self.nuclfile, self.chromfilefwd, rev, outfile))
        self.assertEqual(result, 0)
        with open(outfile, 'r') as resfile:
            arr = resfile.readlines()[0].split(',')
        self.assertEqual(len(arr), 147)
        self.assertEqual(arr[0].strip(), '270.0')
        self.assertEqual(arr[-1].strip(), '297.0')

    def test_fwd1(self):
        outfile = self.configurations['fwd1']['outfile']
        rev = self.configurations['fwd1']['rev']
        result = os.system("python3 {} {} {} {} {}".format(self.scriptname, self.nuclfile, self.chromfilefwd, rev, outfile))
        self.assertEqual(result, 0)
        with open(outfile, 'r') as resfile:
            arr = resfile.readlines()[0].split(',')

        self.assertEqual(len(arr), 147)
        self.assertEqual(arr[0].strip(), '270.0')
        self.assertEqual(arr[-1].strip(), '363.0')

    def test_rev0(self):
        outfile = self.configurations['rev0']['outfile']
        rev = self.configurations['rev0']['rev']
        result = os.system("python3 {} {} {} {} {}".format(self.scriptname, self.nuclfile, self.chromfilerev, rev, outfile))
        self.assertEqual(result, 0)
        with open(outfile, 'r') as resfile:
            arr = resfile.readlines()[0].split(',')

        self.assertEqual(len(arr), 147)
        self.assertEqual(arr[0].strip(), '264.0')
        self.assertEqual(arr[-1].strip(), '318.0')

    def test_rev1(self):
        outfile = self.configurations['rev1']['outfile']
        rev = self.configurations['rev1']['rev']
        result = os.system("python3 {} {} {} {} {}".format(self.scriptname, self.nuclfile, self.chromfilerev, rev, outfile))
        self.assertEqual(result, 0)
        with open(outfile, 'r') as resfile:
            arr = resfile.readlines()[0].split(',')

        self.assertEqual(len(arr), 147)
        self.assertEqual(arr[0].strip(), '259.0')
        self.assertEqual(arr[-1].strip(), '314.0')

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()