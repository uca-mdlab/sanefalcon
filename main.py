import os
import argparse
import re
import logging
import configparser
from merge import merge_all
from prepare_folders import prepare_train_folder
from prepare_samples import prepare_samples
from nucleosome_detector import create_nucl_files


logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s',
                    filename='sanefalcon.log', filemode='w', level=logging.DEBUG)

logger = logging.getLogger("main")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Launch everything',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('conffile', nargs='?', type=str, default="./sanefalcon.conf", help='path of the configuration file')
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.conffile)

    samtools = config['default']['samtools']
    datafolder = config['default']['datafolder']
    trainfolder = config['default']['trainfolder']
    nucleosomefolder = config['default']['nucleosomefolder']
    nucl_file_template = config['default']['nucltemplate']
    bamlist = config['default']['bamlist']

    logger.info("Starting sanefalcon with configuration file {}".format(args.conffile))
    logger.debug("Data folder = {}".format(datafolder))
    logger.debug("Train folder = {}".format(trainfolder))
    logger.debug("Nucl folder = {}".format(nucleosomefolder))

    prepare_samples(datafolder, trainfolder, samtools)
    logger.info("prepare_samples ok")

    prepare_train_folder(bamlist, datafolder, trainfolder)
    logger.info("prepare_folders ok")

    merge_all(trainfolder)
    logger.info("merge_all ok")

    create_nucl_files(trainfolder, nucl_file_template, 'merge')
    logger.info("nucleosome detector merge ok")
    anti_template = nucl_file_template + '_anti'
    create_nucl_files(trainfolder, anti_template, 'anti')
    logger.info("nucleosome detector anti ok")


