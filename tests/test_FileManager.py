import unittest
import os
from manager.file_manager import FileManager
from manager.utils import Utils
import configparser
import shutil


class TestFileUtils(unittest.TestCase):

    def setUp(self):
        conf_file = './tests/data/test.conf'
        self.config = configparser.ConfigParser()
        self.config.read(conf_file)
        self.fm = FileManager(self.config)
        self.fm.check_paths()

    def test_constructor(self):
        self.assertTrue(os.path.isdir(self.fm.trainfolder))
        self.assertTrue(os.path.isdir(self.fm.rspfolder))
        self.assertTrue(os.path.isdir(self.fm.profilefolder))
        bamlist = Utils.readfile(self.fm.config['training']['bamlist'])
        self.assertEqual(len(bamlist), 2)

    def test_prepare_train_folder(self):
        bamlist = Utils.readfile(self.fm.config['training']['bamlist'])
        batchsize = int(self.fm.config['training']['batchsize'])
        batches = self.fm.prepare_train_folder(bamlist, batchsize)
        self.assertEqual(len(batches), 2)
        self.assertTrue('/tmp/testsuite/b' in batches.keys())
        self.assertTrue(all([len(v) == batchsize for k, v in batches.items()]))

    def tearDown(self):
        try:
            shutil.rmtree(self.fm.trainfolder)
        except FileNotFoundError:
            pass


if __name__ == '__main__':
    unittest.main()
