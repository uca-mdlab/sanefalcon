import unittest
from os.path import isabs
from file_manager import FileManager
import configparser


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
        self.assertEqual(len(manips), 1)
        self.assertTrue(isabs(files[0]))

    def test_prepare_train_folder(self):
        self.fm.prepare_train_folder()

if __name__ == '__main__':
    unittest.main()
