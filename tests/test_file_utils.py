import unittest
from file_manager import Utils


class TestFileUtils(unittest.TestCase):
    def test_readfile_int(self):
        f = './tests/data/chromfile.fwd'
        lines = Utils.readfile(f)
        self.assertEqual(len(lines), 71999)
        self.assertEqual(type(lines[0]), int)

    def test_readfile_str(self):
        f = './tests/data/nuclfile'
        lines = Utils.readfile(f)
        self.assertEqual(len(lines), 173875)
        self.assertEqual(type(lines[0]), str)


if __name__ == '__main__':
    unittest.main()
