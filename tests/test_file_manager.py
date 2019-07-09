import unittest
import os
from file_manager import FileManager
import configparser
import shutil

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

    def test_list_files_to_use(self):
        files, manips = self.fm.list_files_to_use()
        self.assertEqual(len(files), 10)
        self.assertEqual(len(manips), 5)
        self.assertTrue(os.path.isabs(files[0]))

    def test_prepare_train_folder(self):
        self.fm.prepare_train_folder()
        self.assertEqual(len(os.listdir(self.fm.trainfolder)), 3)
        links = []
        for root, subdir, files in os.walk(self.fm.trainfolder):
            for f in files:
                links.append(os.path.join(root, f))
        self.assertTrue(all([os.path.islink(x) for x in links]))

    def tearDown(self):
        try:
            shutil.rmtree(self.fm.trainfolder)
        except FileNotFoundError:
            pass


if __name__ == '__main__':
    unittest.main()
