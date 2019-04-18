import os
import string
import re
import logging

logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s',
                    filename='sanefalcon.log', filemode='w', level=logging.DEBUG)

logger = logging.getLogger(__name__)


CHROMOSOMES = range(1, 23)


class Utils:
    @staticmethod
    def readfile(file):
        """
        Read a file and returns a list of stripped lines.
        Tries 'int' conversion for read-start positions
        :param file: the file to read
        :return: a list of either 'str' or 'int'
        """
        with open(file) as f:
            lines = [x.strip() for x in f.readlines()]
        f.close()

        try:
            lines = list(map(int, lines))
        except ValueError:
            pass

        return lines

    @staticmethod
    def prepare_batches(lis, n):
        # Yield successive n-sized chunks from manip_list.
        for i in range(0, len(lis), n):
            yield lis[i:i + n]


class FileManager:
    def __init__(self, config):
        self.datafolder = os.path.abspath(config['default']['datafolder'])
        self.trainfolder = os.path.abspath(config['default']['trainfolder'])
        self.profilefolder = os.path.abspath(config['default']['profilefolder'])
        self.batch_size = int(config['default']['batchsize'])
        self.rspfolder = os.path.abspath(config['default']['rspfolder'])
        self.nucl_file_template = config['default']['nucltemplate']
        self.anti_file_template = self.nucl_file_template + '_anti'
        self.bamlist = Utils.readfile(os.path.abspath(config['default']['bamlist']))
        self.manips = {}
        self.merge_file_lists = {}
        logger.debug("Data folder = {}".format(self.datafolder))
        logger.debug("Train folder = {}".format(self.trainfolder))

    def list_files_to_use(self):
        if not self.bamlist:
            exit('Unable to find the list of bam files to use. Aborting.')
        list_files = [os.path.basename(x) for x in self.bamlist]
        files_to_link = []
        manip_list = set()
        for root, subdir, files in os.walk(self.datafolder):
            for f in files:
                if f in list_files:
                    files_to_link.append(os.path.join(root, f))
                    manip_list.add(os.path.split(root)[1])

        return files_to_link, list(manip_list)

    def prepare_train_folder(self):
        """
        create the train folder by symlinking the original .bam and .bai files into a, b, c, .. subfolders
        """
        letters = list(string.ascii_lowercase)
        files_to_link, manip_list = self.list_files_to_use()
        logger.debug('Found {} files_to_link in {} manips'.format(len(files_to_link), len(manip_list)))
        batches = {}
        for num_batch, batch in enumerate(Utils.prepare_batches(manip_list, self.batch_size)):
            batches[letters[num_batch]] = batch
            logger.debug("Batch {} ({}): {}".format(letters[num_batch], len(batch), batch))

        for batch_name, batch_list in batches.items():
            logger.debug("Preparing Batch {} with {} items".format(batch_name, len(batch_list)))
            workingdir = os.path.join(self.trainfolder, batch_name)
            try:
                os.mkdir(workingdir)
                logger.debug("Created folder: {}".format(workingdir))
            except FileExistsError:
                logger.warning("Folder {} exists, skipping...".format(workingdir))

            for manip in batch_list:
                manip_regex = re.compile(manip)
                files = [f for f in files_to_link if re.search(manip_regex, f)]
                logger.debug("Batch {}: {}".format(batch_name, files))
                for fname in files:
                    run = os.path.split(os.path.split(fname)[0])[1]
                    runpath = os.path.join(workingdir, run)
                    if not os.path.isdir(runpath):
                        os.mkdir(runpath)
                    logger.debug("runpath {} exists ".format(runpath))
                    os.symlink(fname, os.path.join(runpath, os.path.split(fname)[1]))
                    bai_fname = fname + ".bai"
                    os.symlink(bai_fname, os.path.join(runpath, os.path.split(bai_fname)[1]))

        logger.info("Batches created with symlinks")

    def find_all_manips_per_subfolder(self):
        subfolders = [f.path for f in os.scandir(self.trainfolder) if f.is_dir()]
        for subdir in subfolders:
            self.manips[subdir] = [os.path.basename(f.path) for f in os.scandir(subdir) if f.is_dir()]

    def prepare_merge_file_lists(self):
        """
        lists of file for merge.py script
        :return: {'chr1': {'a': [chr1.start.fwd, chr1.start.rev], 'b': [...]}}
        """
        all_start_files = [os.path.join(self.rspfolder, f) for _, _, files in os.walk(self.rspfolder) for f in files]

        files_dic = dict.fromkeys(list(map(str, CHROMOSOMES)), list())
        for chrom in CHROMOSOMES:
            string_pattern = "\.bam\.{}\.start".format(chrom)
            pattern = re.compile(string_pattern)
            files_dic[str(chrom)] = list(filter(lambda x: re.search(pattern, x), all_start_files))

        if not self.manips:
            self.find_all_manips_per_subfolder()
            logger.debug('Loaded manips per subfolder in {}'.format(self.trainfolder))

        for chrom, chrom_start_files in files_dic.items():
            logger.debug('Preparing start_files from chrom {}'.format(chrom))
            files_per_subdir = {}
            for subdir, manip_names in self.manips.items():
                start_per_manip = []
                for n in manip_names:
                    start_per_manip.extend(list(filter(lambda x: re.search(n, x), chrom_start_files)))
                if subdir in files_per_subdir:
                    files_per_subdir[subdir].extend(start_per_manip)
                else:
                    files_per_subdir[subdir] = start_per_manip
            self.merge_file_lists[chrom] = files_per_subdir

    def get_merge_file_lists(self):
        if not self.merge_file_lists:
            self.prepare_merge_file_lists()
        return self.merge_file_lists

    def find_merge_files_in_subdirectories(self):
        merge_files = []
        pattern_name = re.compile("merge.\d{1,2}")
        pattern_subdir = re.compile("/[a-z]/")
        for root, subdirs, files in os.walk(self.trainfolder):
            for fname in files:
                filename = os.path.join(root, fname)
                if re.match(pattern_name, fname) and re.search(pattern_subdir, filename):
                    merge_files.append(filename)
        return merge_files

    def find_merge_anti_files(self):
        merge_files = []
        anti_files = []
        root_merge_files = []
        pattern_name = re.compile("merge.\d{1,2}")
        anti_pattern_name = re.compile("anti.\d{1,2}")
        pattern_subdir = re.compile("/[a-z]/")
        for root, subdirs, files in os.walk(self.trainfolder):
            for fname in files:
                filename = os.path.join(root, fname)
                if re.match(pattern_name, fname) and re.search(pattern_subdir, filename):
                    merge_files.append(filename)
                elif re.match(pattern_name, fname) and not re.search(pattern_subdir, filename):
                    root_merge_files.append(filename)
                if re.match(anti_pattern_name, fname) and re.search(pattern_subdir, filename):
                    anti_files.append(filename)
        return merge_files, anti_files, root_merge_files


if __name__ == '__main__':
    import configparser

    config = configparser.ConfigParser()
    config.read('tests/data/test.conf')
    f = FileManager(config)

    a, b = f.list_files_to_use()
    print(f.bamlist)

