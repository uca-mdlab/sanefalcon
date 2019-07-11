import os
import string
import re
import logging
from prepare_samples import prepare_samples
from collections import defaultdict

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
        self.samtools = config['default']['samtools']
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
        self.rspfiles = defaultdict(dict)
        logger.debug("Data folder = {}".format(self.datafolder))
        logger.debug("Train folder = {}".format(self.trainfolder))
        self.check_conf()

    def check_conf(self):
        if not self.bamlist:
            exit('Unable to find the list of bam files to use. Aborting.')
        for destdir in [self.trainfolder, self.rspfolder, self.profilefolder]:
            if not os.path.isdir(destdir):
                logger.warning('{} not found. Creating...'.format(destdir))
                os.makedirs(destdir)

    def prepare_train_folder(self):
        letters = list(string.ascii_lowercase)
        batches = defaultdict(list)
        logger.debug('Found {} files to link'.format(len(self.bamlist)))
        logger.info('Preparing train folder: symlinking bam files...')
        for num_batch, batch in enumerate([l for l in Utils.prepare_batches(self.bamlist, self.batch_size)]):
            for filename in batch:
                basename = os.path.basename(filename)
                batch_dir = os.path.join(self.trainfolder, letters[num_batch])
                try:
                    os.makedirs(batch_dir)
                    logger.debug("Created folder: {}".format(batch_dir))
                except FileExistsError:
                    logger.warning("Folder {} exists, skipping...".format(batch_dir))
                link_name = os.path.join(batch_dir, basename)
                try:
                    os.symlink(filename, link_name)
                    os.symlink(filename + '.bai', link_name + '.bai')
                except FileExistsError:
                    logger.warning("Link {} exists, skipping...".format(link_name))
                batches[batch_dir].append(link_name)
        logger.info("Batches created with symlinks")

        logger.info('Preparing train folder: prepare_samples...')
        rspfiles = prepare_samples(self.datafolder, self.rspfolder, self.samtools)
        for letter, list_of_fnames in batches.items():
            bam_names = [os.path.basename(f) for f in list_of_fnames]
            for bam in bam_names:
                rsp = list(filter(lambda x: re.search(bam, x), rspfiles))
                self.rspfiles[letter][bam] = rsp

        logger.info('Preparing train folder: mapping complete.')
        return self.rspfiles

    def get_start_files_per_subdir(self):
        merge_file_list = defaultdict(dict)
        for subdir, dic in self.rspfiles.items():
            fwd = defaultdict(list)
            rev = defaultdict(list)
            for bamname, rspfiles in dic.items():
                fwds = [f for f in rspfiles if f.endswith('.fwd')]
                revs = [f for f in rspfiles if f.endswith('.rev')]
                for chrom in range(1, 23):
                    pattern = re.compile(r'\.{}\.'.format(chrom))
                    fwd[chrom].extend(list(filter(lambda x: re.search(pattern, x), fwds)))
                    rev[chrom].extend(list(filter(lambda x: re.search(pattern, x), revs)))

            res = {'fwd': fwd, 'rev': rev}
            merge_file_list[subdir] = res
        return merge_file_list

    def find_all_manips_per_subfolder(self):
        subfolders = [f.path for f in os.scandir(self.trainfolder) if f.is_dir()]
        for subdir in subfolders:
            self.manips[subdir] = [os.path.basename(f.path) for f in os.scandir(subdir) if f.is_dir()]
        for path_to_exclude in [self.profilefolder, self.rspfolder]:  # just in case FIXME
            if path_to_exclude in self.manips:
                del self.manips[path_to_exclude]
        return self.manips

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
    b = f.prepare_train_folder()

