import unittest
import os
from file_manager import FileManager
import configparser
import shutil
import warnings


class TestFileUtils(unittest.TestCase):

    def setUp(self):
        conf_file = './tests/data/test.conf'
        self.config = configparser.ConfigParser()
        self.config.read(conf_file)
        self.fm = FileManager(self.config)

    def test_constructor(self):
        self.assertIsNotNone(self.fm.datafolder)
        self.assertIsNotNone(self.fm.trainfolder)
        self.assertIsNotNone(self.fm.profilefolder)
        self.assertIsInstance(self.fm.batch_size, int)
        self.assertIsNotNone(self.fm.rspfolder)
        self.assertIsNotNone(self.fm.nucl_file_template)
        self.assertEqual(self.fm.nucl_file_template + '_anti', self.fm.anti_file_template)
        self.assertIsNotNone(self.fm.bamlist)
        self.assertIsInstance(self.fm.manips, dict)
        self.assertIsInstance(self.fm.merge_file_lists, dict)

    def test_prepare_train_folder(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = self.fm.prepare_train_folder()
            self.assertTrue(os.path.isdir(self.fm.trainfolder))
            self.assertTrue(os.path.isdir(self.fm.rspfolder))
            self.assertEqual(len(res.keys()), 2)
            self.assertEqual(len(res['a']['p1.bam']), 44)

    def tearDown(self):
        try:
            shutil.rmtree(self.fm.trainfolder)
        except FileNotFoundError:
            pass


if __name__ == '__main__':
    unittest.main()
