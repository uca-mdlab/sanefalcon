import os
import logging
import string
from collections import defaultdict, OrderedDict
from .utils import Utils
from log_setup import setup_logger
import numpy as np


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

    def prepare_train_folder(self, bamlist):
        letters = list(string.ascii_lowercase)
        batches = defaultdict(list)
        logger.debug('Found {} files to link'.format(len(bamlist)))

        reads_count = {}
        for bamfile in bamlist:
            _, num_count = Utils.count_number_of_reads_per_bam(bamfile)
            reads_count[bamfile] = num_count
        tot_reads = sum(reads_count.values())
        logger.debug(f'Total number of reads to balance: {tot_reads}')

        # logger.info('Preparing train folder: symlinking bam files...')

        try:
            for root, dirs, files in os.walk(self.trainfolder):
                for my_file in files:
                    full_path = os.path.join(root, my_file)
                    batches[os.path.dirname(full_path)].append(full_path)
            logger.debug('Batches retrieved:')
            for k, v in batches.items():
                logger.debug(f"{k}, {len([x for x in v if x.endswith('.bam')])}")
        except FileNotFoundError:
            logger.info('Training dir empty, creating from scratch')

        if not batches:
            runs = defaultdict(list)
            for bam in bamlist:
                run = os.path.basename(bam).split('_')[0]
                runs[run].append(bam)

            ordered = OrderedDict(sorted(runs.items(), key=lambda x: chr(int(x[0][1:]))))

            l_ord = list(ordered)
            logger.info(f'Found {len(ordered)} distinct runs')
            result = defaultdict(list)

            low_n, high_n = Utils.compute_num_batches(len(ordered))
            # num_batches = high_n
            num_batches = low_n   # FIXME remove
            logger.info(f'num batches = {num_batches}')

            rows = []
            for b in Utils.prepare_batches(l_ord, num_batches):
                rows.append(list(b))
            res = list(np.array(rows).T)  # transpose to avoid subsequent runs in the same batch
            for i, batch in enumerate(res):
                logger.info(f'Batch {letters[i]}, runs: {batch}, - num samples = '
                            f'{sum(len(ordered[k]) for k in batch)}')
                result[letters[i]].extend(ordered[k] for k in batch)

            result = Utils.balance_batches(result, reads_count)

            for batch_name, l in result.items():
                batch_dir = os.path.join(self.trainfolder, batch_name)
                os.makedirs(batch_dir)
                flatten = [item for sublist in l for item in sublist]
                bam_links = [os.path.join(batch_dir, os.path.basename(item)) for item in flatten]
                bai_links = [os.path.join(batch_dir, os.path.basename(item) + '.bai') for item in flatten]
                for orig, bam_link, bai_link in zip(flatten, bam_links, bai_links):
                    try:
                        os.symlink(orig, bam_link)
                        os.symlink(orig + '.bai', bai_link)
                    except FileExistsError:
                        pass
                    except FileNotFoundError:
                        logger.error(f'{orig} not found??')
                batches[batch_dir] = bam_links

        return batches

    def create_fake_batch_for_testing(self, bamlist):
        batch = [os.path.join(self.datafolder, x) for x in bamlist]
        return {self.testfolder: batch}

