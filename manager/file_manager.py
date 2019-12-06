import os
import logging
import string
from collections import defaultdict
from .utils import Utils
from log_setup import setup_logger

logger = setup_logger(__name__, 'logs/manager.log')


class FileManager:

    def __init__(self, config):
        self.datafolder = config['folders']['data']
        self.trainfolder = config['folders']['train']
        self.profilefolder = config['folders']['profiles']
        self.rspfolder = config['folders']['rsp']
        self.testfolder = config['folders']['test']
        self.testprofilefolder = config['folders']['testprofiles']
        self.config = config

    def check_paths(self):
        for name, path in self.__dict__.items():
            if name in ['datafolder', 'config']:
                continue
            if not os.path.isdir(path):
                logger.warning(f"{name} not found. Creating...")
                os.makedirs(path)

    def get_data(self, mask=''):
        return Utils.get_data(self.datafolder, mask)

    def get_train_data(self, mask=''):
        return Utils.get_data(self.trainfolder, mask)

    def get_profile_data(self, mask=''):
        return Utils.get_data(self.profilefolder, mask)

    def get_read_start_positions_data(self, mask=''):
        return Utils.get_data(self.rspfolder, mask)

    def prepare_train_folder(self, bamlist, batch_size):
        letters = list(string.ascii_lowercase)
        batches = defaultdict(list)
        logger.debug('Found {} files to link'.format(len(bamlist)))
        logger.info('Preparing train folder: symlinking bam files...')
        for num_batch, batch in enumerate([l for l in Utils.prepare_batches(bamlist, batch_size)]):
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
                except FileNotFoundError:
                    pass
                batches[batch_dir].append(link_name)
        logger.info("Batches created with symlinks")
        return batches

    def create_fake_batch_for_testing(self, bamlist):
        batch = [os.path.join(self.datafolder, x) for x in bamlist]
        return {self.testfolder: batch}

