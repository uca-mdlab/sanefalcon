import unittest
from manager.utils import Utils


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

    def test_preparebatches(self):
        l = ['a', 'b', 'c', 'd', 'e', 'f']
        for i in range(1, len(l)):
            num_odd_chunks = 0
            res = list(Utils.prepare_batches(l, i))
            try:
                # evenly sized chunks
                self.assertTrue(all([len(x) == i for x in res]))
            except:
                # only 1 chunk has different length
                num_odd_chunks += 1
            if num_odd_chunks:
                self.assertEqual(num_odd_chunks, 1)


if __name__ == '__main__':
    unittest.main()
