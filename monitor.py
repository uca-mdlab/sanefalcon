import configparser
import logging
import os
import re

from collections import defaultdict

from manager.file_manager import FileManager
from manager.utils import Utils
from manager.rsp import RspBuilder

logging.basicConfig(format='%(asctime)s - %(name)-21s - %(levelname)-8s: %(message)s',
                    filename='sanefalcon.log', filemode='w', level=logging.DEBUG)

logger = logging.getLogger('monitor')


def get_rsp_batches_mapping(batches, rspfiles):
    mapping = defaultdict(dict)
    for subdir, list_of_fnames in batches.items():
        bam_names = [os.path.basename(f) for f in list_of_fnames]
        for bam in bam_names:
            rsp = list(filter(lambda x: re.search(bam, x), rspfiles))
            mapping[subdir][bam] = rsp

    logger.info('Mapping complete.')
    return mapping


if __name__ == '__main__':
    config = configparser.ConfigParser()
    config.read('local.conf')
    logger.info('Starting...')
    f = FileManager(config)
    f.check_paths()
    rs = RspBuilder(config)
    # print(f.get_data())
    # print(f.get_data(mask='.bam$'))
    bamlist = Utils.readfile(config['training']['bamlist'])
    batchsize = int(config['training']['batchsize'])
    batches = f.prepare_train_folder(bamlist, batchsize)
    rspfiles = rs.prepare_samples(f.datafolder, f.rspfolder)
    mapping = get_rsp_batches_mapping(batches, rspfiles)
    print(mapping)