import unittest
from os.path import isabs
from file_manager import FileManager
import configparser


class TestFileUtils(unittest.TestCase):

    def setUp(self):
        conf_file = './tests/data/test.conf'
        self.config = configparser.ConfigParser()
        self.config.read(conf_file)

    def test_constructor(self):
        fm = FileManager(self.config)
        self.assertIsNotNone(fm.datafolder)
        self.assertIsNotNone(fm.trainfolder)
        self.assertIsNotNone(fm.profilefolder)
        self.assertIsInstance(fm.batch_size, int)
        self.assertIsNotNone(fm.rspfolder)
        self.assertIsNotNone(fm.nucl_file_template)
        self.assertEqual(fm.nucl_file_template + '_anti', fm.anti_file_template)
        self.assertIsNotNone(fm.bamlist)
        self.assertIsInstance(fm.manips, dict)
        self.assertIsInstance(fm.merge_file_lists, dict)

    def test_list_files_to_use(self):
        fm = FileManager(self.config)
        files, manips = fm.list_files_to_use()
        self.assertEqual(len(files), 10)
        self.assertEqual(len(manips), 1)
        self.assertTrue(isabs(files[0]))


if __name__ == '__main__':
    unittest.main()
