import os
import argparse
import re
import logging
import configparser
from merge import merge_all
from prepare_folders import prepare_train_folder
from prepare_samples import prepare_samples
from nucleosome_detector import create_nucleosome_files
from getProfileParallel import get_data, submit_process
import multiprocessing as mp


logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s',
                    filename='sanefalcon.log', filemode='w', level=logging.DEBUG)

logger = logging.getLogger("main")


def run_profileParallel(trainfolder, nucl_stub):
    outfolder = os.path.join(trainfolder, "profiles")
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
        logger.info('run_profileParallel: Created out folder {}'.format(outfolder))

    data = get_data(trainfolder, outfolder, nucl_stub)  # all the available data
    input_list = [(chrom, dic, outfolder) for chrom, dic in data.items()]
    print(data,"DATA")
    logger.info("Launching multiprocessing pool...")
    num_cores = mp.cpu_count()
    with mp.Pool(num_cores) as pool:
        finished = pool.map(submit_process, input_list)
        logger.info("Done. Result = {}".format(len(finished) == 22))


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
    nucl_file_template = config['default']['nucltemplate']
    anti_file_template = nucl_file_template + '_anti'
    bamlist = config['default']['bamlist']

    logger.info("Starting sanefalcon with configuration file {}".format(args.conffile))
    logger.debug("Data folder = {}".format(datafolder))
    logger.debug("Train folder = {}".format(trainfolder))

    # prepare_samples(datafolder, trainfolder, samtools)
    # logger.info("prepare_samples ok")
    #
    # try:
    #     prepare_train_folder(bamlist, datafolder, trainfolder)
    #     logger.info("prepare_folders ok")
    # except FileExistsError:
    #     logger.info("train folder symlinks already in place")
    #
    # merge_all(trainfolder)
    # logger.info("merge_all ok")
    #
    # create_nucleosome_files(trainfolder,nucl_file_template,anti_file_template)
    # logger.info("nucleosome ok")
    print(trainfolder,"TRAINFOLDER ++++++++++")
    run_profileParallel(trainfolder, nucl_file_template)
    logger.info("run profile parallel ok")
